"""Explicit mass coordinates and dual-peak uncertainty; independent of candidate ranking."""
from dataclasses import dataclass
import math
from knowledge_source import mass_constants

@dataclass(frozen=True)
class MassCoordinate:
    kind: str
    value: float
    tolerance_da: float
    mass_basis: str
    charge: int = 0
    intrinsic_charge: int = 0
    ions: tuple = ()
    isotope_shift: float = 0.0
    assumption: str = 'explicit'
    tolerance_ppm: float | None = None

    @property
    def multiplier(self):return abs(self.charge) if self.kind=='mz' else 1

    @property
    def offset(self):
        c=mass_constants();e=c['electron_da'];p=c['proton_da']
        atom_table=c['monoisotopic_atoms' if self.mass_basis=='monoisotopic' else 'average_atoms']
        ions={'H':p,'-H':-p,'Na':atom_table['Na']-e,'K':atom_table['K']-e}
        return sum(ions[name]*count for name,count in self.ions)-self.intrinsic_charge*e+self.isotope_shift

    def predict(self,element_mass):return (element_mass+self.offset)/self.multiplier
    def elemental_reference(self):return self.value*self.multiplier-self.offset
    def elemental_tolerance(self):return self.tolerance_da*self.multiplier
    def match(self,element_mass):
        predicted=self.predict(element_mass);error=self.value-predicted
        ppm=1e6*error/predicted if predicted>0 else None
        match=abs(error)<=self.tolerance_da if self.tolerance_ppm is None else ppm is not None and abs(ppm)<=self.tolerance_ppm
        return {'match':match,'predicted':predicted,'error_da':error,
                'error_ppm':ppm,
                'assumption':self.assumption,'mass_basis':self.mass_basis}

def finite(value,label,positive=False):
    if type(value) not in (float,int) or not math.isfinite(value) or positive and value<=0:
        raise ValueError(label+' must be a finite '+('positive number' if positive else 'number'))
    return value

def parse_coordinate(raw,assumption='explicit'):
    allowed={'kind','value','tolerance','tolerance_unit','mass_basis','charge','intrinsic_charge','ions','isotope_shift'}
    if not isinstance(raw,dict) or set(raw)-allowed:raise ValueError('Unknown observation fields')
    kind=raw.get('kind');value=finite(raw.get('value'),'observation value',True)
    tolerance=finite(raw.get('tolerance'),'tolerance',True);unit=raw.get('tolerance_unit')
    if unit not in {'Da','ppm'}:raise ValueError('Tolerance unit must be Da or ppm')
    tol=tolerance if unit=='Da' else value*tolerance/1e6
    if kind not in {'neutral_mono','neutral_average','mz'}:raise ValueError('Unsupported mass coordinate')
    basis=raw.get('mass_basis','average' if kind=='neutral_average' else 'monoisotopic')
    if basis not in {'monoisotopic','average'}:raise ValueError('Explicit monoisotopic or average basis required')
    if kind=='neutral_mono' and basis!='monoisotopic' or kind=='neutral_average' and basis!='average':
        raise ValueError('Coordinate kind and mass basis disagree')
    isotope=finite(raw.get('isotope_shift',0),'isotope shift')
    if isotope<0 or isotope and basis=='average':raise ValueError('Resolved isotope assignment requires a nonnegative monoisotopic shift')
    if kind!='mz':
        if {'charge','intrinsic_charge','ions'}&raw.keys():raise ValueError('Neutral mass must not carry ion components')
        return MassCoordinate(kind,value,tol,basis,isotope_shift=isotope,assumption=assumption,tolerance_ppm=tolerance if unit=='ppm' else None)
    q=raw.get('charge');intrinsic=raw.get('intrinsic_charge',0);ions=raw.get('ions')
    if type(q) is not int or q==0 or type(intrinsic) is not int:raise ValueError('m/z needs explicit nonzero integer charge and integer intrinsic charge')
    if not isinstance(ions,dict) or any(k not in {'H','-H','Na','K'} or type(v) is not int or v<=0 for k,v in ions.items()):
        raise ValueError('Explicit H/-H/Na/K ion counts required')
    if 'H' in ions and '-H' in ions:raise ValueError('Cannot add and remove protons simultaneously')
    charge=intrinsic+sum((-1 if name=='-H' else 1)*n for name,n in ions.items())
    if q!=charge:raise ValueError('Ion plus intrinsic charge does not equal observed charge')
    return MassCoordinate(kind,value,tol,basis,q,intrinsic,tuple(sorted(ions.items())),isotope,assumption,tolerance if unit=='ppm' else None)

def hypotheses(raw):
    """Finite explicit alternatives. No undeclared isotope, adduct, or charge guessing."""
    if 'ion_hypotheses' not in raw:return [parse_coordinate(raw)]
    options=raw['ion_hypotheses']
    if not isinstance(options,list) or not options or len(options)>32:raise ValueError('Declare 1–32 finite ion hypotheses')
    base={k:v for k,v in raw.items() if k!='ion_hypotheses'}
    if {'charge','ions','intrinsic_charge'}&base.keys():raise ValueError('Declare ions either on observation or in hypotheses')
    out=[]
    for i,option in enumerate(options):
        if not isinstance(option,dict) or set(option)-{'charge','ions','intrinsic_charge','isotope_shift'}:
            raise ValueError('Unknown ion hypothesis fields')
        out.append(parse_coordinate({**base,**option},assumption=f'declared_hypothesis_{i+1}'))
    return list(dict.fromkeys(out))

def formula_mass(formula,basis='monoisotopic'):
    table=mass_constants()['monoisotopic_atoms' if basis=='monoisotopic' else 'average_atoms']
    if basis not in {'monoisotopic','average'}:raise ValueError('Unsupported mass basis')
    if not isinstance(formula,dict) or any(e not in table or type(n) is not int for e,n in formula.items()):
        raise ValueError('Unsupported element or non-integer count')
    return sum(table[e]*n for e,n in formula.items())

def peak_pair(first,second,error_model='bounded'):
    if first.mass_basis!=second.mass_basis:raise ValueError('Peak mass bases disagree')
    a,b=first.elemental_tolerance(),second.elemental_tolerance()
    if error_model=='bounded':error=a+b
    elif error_model=='independent_gaussian':error=math.hypot(a,b)
    else:raise ValueError('Use bounded errors or explicitly independent Gaussian standard uncertainties')
    return {'delta_element_mass':second.elemental_reference()-first.elemental_reference(),
            'tolerance_da':error,'mass_basis':first.mass_basis,'error_model':error_model,
            'first_assumption':first.assumption,'second_assumption':second.assumption,
            'reference_kind':'neutral_element_sum' if not first.intrinsic_charge and not second.intrinsic_charge else 'declared_ion_element_reference'}

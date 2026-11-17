"""Conservative Chinese/English peptide input parsing; missing precision or structure stays explicit."""
import re
from case_input import read_sequence

ALIASES={'哌啶':'piperidine','piperidine':'piperidine','二甲亚砜':'dmso','dmso':'dmso',
         'edt':'edt','tfa':'acid_exposure','甲醛':'formaldehyde','甲酸':'formic_acid','氧化':'oxidation',
         '三苯基膦':'triphenylphosphine','锌':'zinc','碱处理':'base_exposure','存放':'storage'}
ALIASES.update({'三氟乙酸':'tfa','tfa':'tfa','甲醇':'methanol','methanol':'methanol','meoh':'methanol',
    '乙酸':'acetic_acid','acetic acid':'acetic_acid','formic acid':'formic_acid','二甲基亚砜':'dmso',
    '二异丙基乙胺':'diea','dipea':'diea','diea':'diea','dcc':'dcc','dic':'dic','hobt':'hobt',
    'bop':'bop','hbtu':'hbtu','mtbe':'mtbe','甲基叔丁基醚':'mtbe','dtt':'dtt','tcep':'tcep',
    '二硫苏糖醇':'dtt','乙二硫醇':'edt','三异丙基硅烷':'tis','tis':'tis','dmf':'dmf',
    'dcm':'dcm','二氯甲烷':'dcm','乙腈':'acn','acn':'acn','肼':'hydrazine',
    '羟胺':'hydroxylamine','甲醛':'formaldehyde','丙酮':'acetone'})
NUMBER=r'[+−-]?\d+(?:\.\d+)?'
def parse_text(text):
    if not isinstance(text,str) or not text.strip():raise ValueError('Provide the original input text')
    missing=[];seq=re.search(r'(?:序列|sequence)\s*[:：为是]?\s*([A-Za-z](?:[A-Za-z]|[\s-](?=[A-Za-z]))*)',text,re.I)
    if not seq:missing.append('多肽序列');sequence=None
    else:
        sequence=seq[1].strip().replace(' ','').replace('-','')
        try:read_sequence(sequence)
        except ValueError:missing.append('无歧义的一字母序列');sequence=None
    nterm='free' if re.search(r'N\s*端\s*(?:游离|free|未保护)',text,re.I) else 'acetyl' if re.search(r'N\s*端\s*(?:乙酰|Ac)',text,re.I) else None
    cterm='amide' if re.search(r'C\s*端\s*(?:酰胺|amide|NH2)',text,re.I) else 'acid' if re.search(r'C\s*端\s*(?:游离酸|羧酸|酸|acid|COOH)',text,re.I) else None
    if nterm is None:missing.append('N端状态')
    if cterm is None:missing.append('C端状态')
    tol=re.search(r'(?:容差|误差|tolerance)\s*[:：±=]?\s*('+NUMBER+r')\s*(Da|ppm)',text,re.I)
    tolerance=float(tol[1].replace('−','-')) if tol else None;unit=tol[2] if tol else None
    if tolerance is None or tolerance<=0:missing.append('实际 Da／ppm 容差')
    unit='Da' if unit and unit.lower()=='da' else unit.lower() if unit else None
    delta=re.search(r'(?:Δ\s*M|delta\s*M|质量差)\s*[:：=]?\s*('+NUMBER+r')\s*Da',text,re.I)
    amount=re.search(r'(?:比理论|相对理论)\s*(多(?:出)?|少|减少|增加)\s*('+NUMBER+r')\s*Da',text,re.I)
    mz=re.search(r'm\s*/\s*z\s*[:：=]?\s*('+NUMBER+r')',text,re.I)
    neutral=re.search(r'(?:中性(?:单同位素)?质量|neutral\s+mass)\s*[:：=为]?\s*('+NUMBER+r')',text,re.I)
    observation=None
    if delta or amount:
        value=float((delta[1] if delta else amount[2]).replace('−','-'))
        if amount and amount[1] in {'少','减少'}:value=-abs(value)
        if unit and unit!='Da':missing.append('质量差的 Da 误差；不能把相对差值 ppm 当成峰质量误差')
        observation={'kind':'delta_approx' if re.search('低分辨|约|approx',text,re.I) else 'delta_mono',
                     'value':value,'reference':'intended_product','tolerance':tolerance,'tolerance_unit':unit}
    elif mz:
        charge=re.search(r'(?:电荷|charge)\s*[:：=为]?\s*([+−-]?\d+)',text,re.I)
        if not charge:
            without_mz=text[:mz.start()]+text[mz.end():]
            charge=re.search(r'\bz\s*[:：=为]?\s*([+−-]?\d+)',without_mz,re.I)
        adduct='Na' if re.search(r'\[M\+Na\]|钠加合',text,re.I) else 'K' if re.search(r'\[M\+K\]|钾加合',text,re.I) else 'H' if re.search(r'质子化|\[M\+(?:\d+)?H\]',text,re.I) else '-H' if re.search(r'去质子化|\[M-(?:\d+)?H\]',text,re.I) else None
        q=int(charge[1].replace('−','-')) if charge else None
        if q is None or q==0:missing.append('m/z 的非零整数电荷')
        if adduct is None:missing.append('m/z 的加合离子，或允许枚举的有限假设')
        if q and adduct:
            count=abs(q) if adduct in {'H','-H'} else 1
            observation={'kind':'mz','value':float(mz[1]),'charge':q,'ions':{adduct:count},'tolerance':tolerance,'tolerance_unit':unit}
    elif neutral:
        basis='average' if re.search('平均质量|average',text,re.I) else 'monoisotopic' if re.search('单同位素|monoisotopic|HRMS|高分辨',text,re.I) else None
        if basis is None:missing.append('中性质量的平均／单同位素口径')
        observation={'kind':'neutral_average' if basis=='average' else 'neutral_mono','value':float(neutral[1]),'tolerance':tolerance,'tolerance_unit':unit}
    else:missing.append('有明确标签的质量差、中性质量或 m/z')
    present=[];absent=[];lower=text.lower()
    for alias,tag in ALIASES.items():
        pattern=(r'(?<![a-z])'+re.escape(alias)+r'(?![a-z])') if alias.isascii() else re.escape(alias)
        for found in re.finditer(pattern,lower):
            before=lower[max(0,found.start()-12):found.start()]
            if re.search(r'(?:未用|没用|没有用|未加入|未加|未使用|没有|没有加|不加|不含|不使用|without|no)\s*(?:过|加入|使用)?\s*$',before):absent.append(tag)
            else:present.append(tag)
    if set(present)&set(absent):missing.append('互相矛盾的试剂信息需澄清')
    if 'tfa' in present:present.append('acid_exposure')
    if 'tfa' in absent:absent.append('acid_exposure')
    if len(list(re.finditer(r'm\s*/\s*z\s*[:：=]?\s*'+NUMBER,text,re.I)))>1 or (mz and (delta or amount or neutral)):
        missing.append('多个质量观测需分别登记；不静默忽略第二个峰')
    if re.search(r'保护|修饰|Fmoc|Boc|Acm|Pbf|Cam|二硫|disulfide|磷酸化',text,re.I):
        missing.append('保护／修饰／连接的明确位点；当前自然语言简式尚不解析该类结构')
    if missing:return {'status':'needs_input','missing':sorted(set(missing)),
                       'parsed':{'sequence':sequence,'nterm':nterm,'cterm':cterm,'observation':observation,
                                 'conditions':{'stage_tags':sorted(set(present)),'stage_tags_absent':sorted(set(absent))}},'input_text':text}
    return {'status':'ready','case':{'case_id':'natural-language','input_text':text,
        'intended_product':{'sequence':sequence,'nterm':nterm,'cterm':cterm},
        'reaction':{'description':text,'stage_tags':sorted(set(present)),'stage_tags_absent':sorted(set(absent))},
        'observations':[observation]}}

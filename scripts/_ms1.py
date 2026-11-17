"""Centroid MS1 TIC/BPC/XIC using indexed mzML; overlap and arbitrary intensity units are explicit."""
import math
from _experimental import number
from case_input import check_keys

def extract_scans(scans,specification):
    check_keys(specification,{'targets','rt_min_seconds','rt_max_seconds'},'MS1 targets')
    targets=specification.get('targets',[])
    if not isinstance(targets,list) or not targets:raise ValueError('Provide one or more explicit target m/z values')
    prepared=[];ids=set()
    for t in targets:
        check_keys(t,{'id','mz','tolerance','tolerance_unit'},'XIC target')
        key=t.get('id');mz=number(t.get('mz'),'target m/z',True);tol=number(t.get('tolerance'),'target tolerance',True)
        if not isinstance(key,str) or not key or key in ids:raise ValueError('Target IDs must be nonempty and unique')
        ids.add(key);unit=t.get('tolerance_unit')
        if unit not in {'Da','ppm'}:raise ValueError('Target tolerance unit must be Da or ppm')
        prepared.append((key,mz,tol if unit=='Da' else mz*tol/1e6))
    low=specification.get('rt_min_seconds',0);high=specification.get('rt_max_seconds')
    number(low,'minimum RT')
    if low<0 or high is not None and number(high,'maximum RT')<low:raise ValueError('RT range must be ordered and nonnegative')
    rows=[];overlaps=0;last=None
    for scan in scans:
        if scan['ms_level']!=1:continue
        if scan.get('centroided') is not True:raise ValueError('MS1 extraction requires explicitly centroided scans; profile data must be centroided first')
        rt=number(scan['rt_seconds'],'scan RT')
        if rt<0 or last is not None and rt<last:raise ValueError('MS1 retention times must be nonnegative and nondecreasing')
        last=rt
        if rt<low or high is not None and rt>high:continue
        peaks=scan.get('peaks',[])
        values={key:0.0 for key,_,_ in prepared};tic=0.0;bpc=0.0
        for peak in peaks:
            if not isinstance(peak,(list,tuple)) or len(peak)!=2:raise ValueError('Each centroid peak needs m/z and intensity')
            mz=number(peak[0],'peak m/z',True);intensity=number(peak[1],'peak intensity')
            if intensity<0:raise ValueError('Centroid intensities must be nonnegative')
            tic+=intensity;bpc=max(bpc,intensity);hits=[]
            for key,target,window in prepared:
                if abs(mz-target)<=window:values[key]+=intensity;hits.append(key)
            if len(hits)>1:overlaps+=1
        rows.append({'rt_seconds':rt,'TIC':tic,'BPC':bpc,'XIC':values})
    auc={key:math.fsum((a['XIC'][key]+b['XIC'][key])*(b['rt_seconds']-a['rt_seconds'])/2 for a,b in zip(rows,rows[1:])) for key,_,_ in prepared}
    return {'module':'ms1_basics','ms1_rows':rows,'target_specification':specification,'area_intensity_seconds':auc,
        'overlapping_peak_target_assignments':overlaps,'no_ms1_in_selected_range':not rows,
        'scope':'raw centroid-window sums and trapezoidal areas; no baseline correction, calibrated concentration or structural assignment'}

def mzml_scans(path):
    import pyopenms as ms
    run=ms.OnDiscMSExperiment()
    if not run.openFile(str(path)):raise ValueError('Indexed mzML required for bounded per-spectrum reads; convert/index the input first')
    for i in range(run.getNrSpectra()):
        s=run.getSpectrum(i)
        if s.getMSLevel()!=1:continue
        mz,intensities=s.get_peaks()
        yield {'ms_level':1,'rt_seconds':s.getRT(),'centroided':s.getType()==ms.SpectrumSettings.SpectrumType.CENTROID,
            'peaks':[(float(x),float(y)) for x,y in zip(mz,intensities)]}

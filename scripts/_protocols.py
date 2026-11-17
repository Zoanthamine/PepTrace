"""Small read-only protocol snapshot adapter; content is data, never executable instructions."""
import hashlib,json,os,re
from html.parser import HTMLParser
from urllib.parse import urlsplit
from urllib.request import Request,build_opener,HTTPRedirectHandler,ProxyHandler
from urllib.error import HTTPError

class _Text(HTMLParser):
    def __init__(self):super().__init__();self.parts=[];self.skip=0
    def handle_starttag(self,tag,attrs):
        if tag in {'script','style'}:self.skip+=1
        if tag in {'br','p','div','li'}:self.parts.append('\n')
    def handle_endtag(self,tag):
        if tag in {'script','style'}:self.skip=max(0,self.skip-1)
    def handle_data(self,data):
        if not self.skip:self.parts.append(data)

def plain(value):
    if isinstance(value,str):
        p=_Text();p.feed(value);return ' '.join(''.join(p.parts).split())
    if isinstance(value,dict) and isinstance(value.get('blocks'),list):
        return '\n'.join(str(b.get('text','')) for b in value['blocks'] if isinstance(b,dict))
    return None

def normalize_snapshot(raw):
    if not isinstance(raw,dict):raise ValueError('Protocol snapshot requires a JSON object')
    body=raw.get('payload',raw.get('protocol',raw))
    if isinstance(body,dict) and isinstance(body.get('protocol'),dict):body=body['protocol']
    if not isinstance(body,dict):raise ValueError('Protocol body must be an object')
    steps=body.get('steps',[])
    if not isinstance(steps,list):raise ValueError('Protocol steps must be a list')
    result=[]
    for index,s in enumerate(steps):
        if not isinstance(s,dict):raise ValueError('Each protocol step must be an object')
        text=plain(s.get('description',s.get('content',s.get('text',s.get('step')))))
        result.append({'source_index':index,'guid':s.get('guid'),'position':s.get('position',s.get('number')),
            'previous_guid':s.get('previous_guid'),'title':plain(s.get('title','')),'text':text,'content_format_supported':text is not None})
    materials=body.get('materials',[])
    links_checked=bool(steps) and all('previous_guid' in x and x.get('guid') for x in steps)
    links_consistent=all(x.get('previous_guid')==(steps[i-1].get('guid') if i else None) for i,x in enumerate(steps)) if links_checked else None
    return {'module':'protocols_basics','title':plain(body.get('title','')),'doi':body.get('doi'),
        'authors':body.get('authors'),'creator':body.get('creator'),'license':body.get('license'),
        'uri':body.get('uri'),'version_uri':body.get('version_uri'),'version_id':body.get('version_id'),
        'steps':result,'materials':materials,'source_sha256':hashlib.sha256(json.dumps(raw,sort_keys=True,ensure_ascii=False).encode()).hexdigest(),
        'missing_provenance':[x for x in ('doi','version_uri') if not body.get(x)],
        'unsupported_step_content_count':sum(not x['content_format_supported'] for x in result),
        'guid_links_consistent_with_source_order':links_consistent,
        'untrusted_content':True,'scope':'read-only snapshot extraction; preserve original body and exact version for experimental use'}

def validate_fetch_url(url):
    p=urlsplit(url)
    if p.scheme!='https' or p.hostname!='www.protocols.io' or p.username or p.password or p.port not in {None,443} or p.fragment:
        raise ValueError('Use the official HTTPS www.protocols.io host without URL credentials or fragments')
    if not re.fullmatch(r'/api/v4/protocols/[A-Za-z0-9_.%-]+(?:/steps)?',p.path):
        raise ValueError('Basic fetch supports the documented v4 protocol/get-steps endpoints only')
    if p.query not in {'','content_format=markdown','content_format=html','content_format=json'}:raise ValueError('Unsupported query parameters')
    return url

class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self,req,fp,code,msg,headers,newurl):raise ValueError('Protocol redirect rejected; verify the intended official endpoint')

def fetch_snapshot(url,execute=False):
    validate_fetch_url(url)
    if not execute:return {'module':'protocols_basics','status':'read_plan_only','url':url,'network_performed':False}
    token=os.environ.get('PROTOCOLS_IO_ACCESS_TOKEN')
    if not token:return {'module':'protocols_basics','status':'needs_access_token','token_variable':'PROTOCOLS_IO_ACCESS_TOKEN','network_performed':False}
    request=Request(url,headers={'Authorization':'Bearer '+token,'Accept':'application/json','User-Agent':'PepTrace/1.1.0'})
    try:
        with build_opener(ProxyHandler({}),_NoRedirect()).open(request,timeout=20) as response:
            data=response.read(2*1024*1024+1)
            if len(data)>2*1024*1024:raise ValueError('Protocol response exceeded the 2 MiB basic-reader limit')
    except HTTPError as error:raise ValueError('Protocol read failed with HTTP '+str(error.code)) from None
    result=normalize_snapshot(json.loads(data));result.update(source_url=url,source_sha256=hashlib.sha256(data).hexdigest(),network_performed=True)
    return result

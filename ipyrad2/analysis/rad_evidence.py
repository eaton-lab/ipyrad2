"""Bounded access to native all-RAD-site evidence in HDF5 format 4.0."""
import json
import h5py
import numpy as np
from ..utils.exceptions import IPyradError
from .evidence import _names


def _selection(indices):
    """Use one HDF5 slice for a sorted contiguous selection; preserve sparse reads."""
    if len(indices) and int(indices[-1])-int(indices[0])+1 == len(indices):
        return slice(int(indices[0]), int(indices[-1])+1)
    return indices


class RadSiteEvidenceReader:
    """Read bases and their zero, one or several native pileup records.

    AD/PL are padded only in the returned selection. `allele_count` and
    `PL_count` distinguish missing native entries from padding. Symbolic allele
    tokens are returned literally; the SNP popgen model must not consume these
    records as if their unseen-allele category were a particular nucleotide.
    """
    def __init__(self,path):
        self.io=h5py.File(path,'r',rdcc_nbytes=128*1024**2)
        try:
            if self.io.attrs.get('format_version')!='4.0' or self.io.attrs.get('evidence_scope')!='rad_sites':
                raise IPyradError('All-locus evidence requires HDF5 format 4.0; assemble with --evidence loci.')
            self.group=self.io['rad_evidence'];g=self.group
            if g.attrs.get('schema_version')!='4.0' or not g.attrs.get('complete',False):
                raise IPyradError('All-RAD evidence is incomplete or unsupported.')
            self.names=_names(g['samples'][:]);self.scaffolds=_names(g['scaffolds'][:]);self.tokens=_names(g['allele_tokens'][:])
            self.status_codes={int(k):v for k,v in json.loads(g.attrs['status_codes']).items()}
            self.nsites=int(g.attrs['sites']);self.nrecords=int(g.attrs['records'])
            ns=len(self.names);nr=self.nrecords
            if g['site_map'].shape!=(self.nsites,3) or g['hard_mask'].shape!=(ns,self.nsites) or g['record_span'].shape!=(self.nsites,2):
                raise IPyradError('All-RAD evidence site shapes are inconsistent.')
            for key in ['DP','GQ','status']:
                if g[key].shape!=(ns,nr):raise IPyradError('All-RAD evidence record shapes are inconsistent.')
            for key in ['allele_offsets','PL_offsets']:
                if g[key].shape!=(nr+1,):raise IPyradError('All-RAD vector offsets are inconsistent.')
            if g['AD'].shape!=(ns,int(g['allele_offsets'][-1])) or g['PL'].shape!=(ns,int(g['PL_offsets'][-1])):
                raise IPyradError('All-RAD vector shapes are inconsistent.')
            if g['allele_ids'].shape!=(int(g['allele_offsets'][-1]),):raise IPyradError('All-RAD allele axes are inconsistent.')
        except BaseException:
            self.io.close();raise

    def __enter__(self):return self
    def __exit__(self,*args):self.io.close()

    def read(self,sites,samples=None):
        raw=np.asarray(sites);sites=np.asarray(sites,dtype=np.int64)
        if sites.ndim!=1 or not np.array_equal(raw,sites) or len(set(sites))!=len(sites) or np.any(sites<0) or np.any(sites>=self.nsites):
            raise IPyradError('RAD site indices must be unique in-range integers.')
        samples=self.names if samples is None else list(samples)
        if len(set(samples))!=len(samples) or any(s not in self.names for s in samples):
            raise IPyradError('RAD evidence requires unique observed sample names.')
        ix=[self.names.index(s) for s in samples];g=self.group
        order=np.argsort(sites);inverse=np.argsort(order);sorted_sites=sites[order]
        spans=g['record_span'][_selection(sorted_sites)][inverse]
        records=np.concatenate([np.arange(int(a),int(a+b),dtype=np.int64) for a,b in spans]) if len(spans) else np.array([],dtype=np.int64)
        ro=np.argsort(records);ri=np.argsort(ro);sr=records[ro]
        d=dict(site_indices=sites,site_map=g['site_map'][_selection(sorted_sites)][inverse],reference=g['reference'][_selection(sorted_sites)][inverse],
               site_status=g['site_status'][_selection(sorted_sites)][inverse],hard_mask=g['hard_mask'][:,_selection(sorted_sites)][ix][:,inverse],
               record_indices=records,record_site=g['record_site'][_selection(sr)][ri],QUAL=g['QUAL'][_selection(sr)][ri])
        for key in ['DP','GQ','status']:d[key]=g[key][:,_selection(sr)][ix][:,ri]
        offsets={}
        for key in ['allele_offsets','PL_offsets']:
            select=np.unique(np.r_[sr,sr+1]);values=g[key][_selection(select)];lookup=dict(zip(select,values))
            offsets[key]=np.array([(lookup[r],lookup[r+1]) for r in records],dtype=np.int64).reshape(-1,2)
        d['allele_count']=np.diff(offsets['allele_offsets'],axis=1).ravel()
        d['PL_count']=np.diff(offsets['PL_offsets'],axis=1).ravel()
        d['alleles']=[]
        for key,okey in [('AD','allele_offsets'),('PL','PL_offsets')]:
            spans=offsets[okey];columns=np.unique(np.concatenate([np.arange(a,b) for a,b in spans])) if len(spans) else np.array([],dtype=np.int64)
            values=g[key][:,_selection(columns)][ix];lookup={int(c):j for j,c in enumerate(columns)}
            counts=np.diff(spans,axis=1).ravel();d[key]=np.full((len(samples),len(records),int(counts.max(initial=0))),-1,dtype='i4')
            ids=g['allele_ids'][_selection(columns)] if key=='AD' else None
            for j,(a,b) in enumerate(spans):
                cols=[lookup[c] for c in range(a,b)];d[key][:,j,:b-a]=values[:,cols]
                if key=='AD':d['alleles'].append(tuple(self.tokens[int(k)] for k in ids[cols]))
        expected=d['allele_count']*(d['allele_count']+1)//2
        if not np.array_equal(expected,d['PL_count']):raise IPyradError('RAD likelihood cardinality disagrees with allele count.')
        active=np.arange(d['PL'].shape[-1])[None,:]<d['PL_count'][:,None]
        d['usable_PL']=(d['status']==0)&(d['DP']>0)&np.all((d['PL']>=0)|~active[None,:,:],axis=-1)
        smap={int(s):i for i,s in enumerate(sites)}
        cell_mask=d['hard_mask'][:,[smap[int(s)] for s in d['record_site']]].astype(bool)
        d['usable_PL'] &= ~cell_mask
        d['is_base_alleles_only']=np.array([all(len(a)==1 and a in 'ACGT' for a in aa) for aa in d['alleles']],dtype=bool)
        return d

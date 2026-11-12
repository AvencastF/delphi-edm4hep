"""Read every EDM4hep frame and report content, identities and relation integrity."""
import collections,csv,json,math,struct,sys
from pathlib import Path
from podio import root_io
rootfile,fadgen,outdir=map(Path,sys.argv[1:]); outdir.mkdir(exist_ok=True)
fad=[]
with fadgen.open('rb') as f:
 while True:
  h=f.read(4); s,=struct.unpack('<i',h); b=f.read(s); assert f.read(4)==h
  n,=struct.unpack_from('<i',b)
  if not n: break
  fad.append([struct.unpack_from('<5i10f',b,4+60*i) for i in range(n)])
reader=root_io.Reader(str(rootfile)); sizes={}; ids=[]; errors=[]; eventrows=[]; truthstats=collections.Counter(); ecs=collections.Counter(); tags=collections.Counter(); max_pdiff=0.; max_vdiff=0.; comparable=0
from prepare_weights import prepare
prepare(Path.cwd())  # independently revalidate the production inputs
with Path('event-weights.csv').open() as stream:
 weight_rows={int(r['hepmc_event']):r for r in csv.DictReader(stream)}
expected_run=-int(json.loads(Path('task.json').read_text())['run'])
observed_weights=[]
reco_diagnostics=dict(no_own_row_charged=0,no_own_row_selected_charged=0,recovered_neutrals=0,charge_codes=collections.Counter(),vecp_charges=collections.Counter(),mass_codes=collections.Counter(),perigee_refitted=collections.Counter(),selection_bits=collections.Counter())
for i,fr in enumerate(reader.get('events')):
 run=int(fr.get_parameter('sDST_EVT_runNumber')); evt=int(fr.get_parameter('sDST_EVT_eventNumber')); ids.append((run,evt))
 weight=float(fr.get_parameter('mc_gen_weight'))
 trial=int(fr.get_parameter('mc_generator_trial'))
 if run!=expected_run or evt not in weight_rows: raise ValueError('Wrong embedded weight identity')
 row=weight_rows[evt]
 if weight!=float(row['weight']) or trial!=int(row.get('generator_trial',-1)): raise ValueError('Embedded weight mismatch')
 observed_weights.append(weight)
 ecs[str(fr.get_parameter('sDST_EVT_ECMS'))]+=1; tags[str(fr.get_parameter('sDST_EVT_dstProcessingTag'))]+=1
 names=set(map(str,fr.getAvailableCollections()))
 for name in names:
  c=fr.get(name); n=c.size()
  d=sizes.setdefault(name,dict(type=str(c.getTypeName()),total=0,nonempty_events=0,min=n,max=n))
  d['total']+=n; d['nonempty_events']+=bool(n); d['min']=min(d['min'],n); d['max']=max(d['max'],n)
 ps=fr.get('sDST_MAIN_Particles'); tr=fr.get('sDST_TRAC_Tracks'); mc=fr.get('sDST_LUJ_GenParticles'); links=fr.get('sDST_TBL_RecoToGen')
 flags=fr.get('sDST_VECP_Particles_SelectionFlag'); status=fr.get('sDST_MAIN_Particles_VecpStatus')
 codes=fr.get('sDST_MAIN_Particles_ChargeCode'); vecq=fr.get('sDST_VECP_Particles_Charge'); masscodes=fr.get('sDST_MAIN_Particles_MassCode')
 for index,particle in enumerate(ps):
  # VecpStatus 1/3: charged PA without its own VECP row, NaN 4-vector; 2: SKELANA's recovered neutral.
  if status[index] in (1,3):
   reco_diagnostics['no_own_row_charged']+=1
   reco_diagnostics['no_own_row_selected_charged']+=int(flags[index]==0)
  reco_diagnostics['recovered_neutrals']+=int(status[index]==2)
  # Raw PA.MAIN charge code and SKELANA's VECP(7) (NaN without an own row), counted as published.
  reco_diagnostics['charge_codes'][int(codes[index])]+=1; reco_diagnostics['vecp_charges'][str(vecq[index])]+=1
  # Raw PA.MAIN mass code; a recovered neutral repeats its source's, so it is counted once.
  if status[index]!=2: reco_diagnostics['mass_codes'][int(masscodes[index])]+=1
 # SKELANA perigee words: 1 refitted (differ from PA.TRAC), 0 equal, -1 no own VECP row.
 for flag in fr.get('sDST_QTRAC_Tracks_PerigeeRefitted'): reco_diagnostics['perigee_refitted'][int(flag)]+=1
 # SKELANA's LVSELE word: which bits (1 = least significant) are set, by VecpStatus; -1 = no own row.
 for word,s in zip(fr.get('sDST_VECP_Particles_SelectionWord'),status):
  word=int(word)&0xFFFFFFFF
  for bit in ([0] if word==0xFFFFFFFF else [b+1 for b in range(32) if word>>b&1]): reco_diagnostics['selection_bits']['status%d_bit%d'%(s,bit)]+=1
 charged=sum(p.getCharge()!=0 for p in ps)
 for index,p in enumerate(ps):
  v=p.getMomentum()
  if all(math.isfinite(x) for x in (v.x,v.y,v.z,p.getEnergy())) == (status[index] in (1,3)): errors.append(['reco_nonfinite' if status[index] not in (1,3) else 'reco_missing_not_nan',evt])
  for t in p.getTracks():
   if not 0<=t.getObjectID().index<tr.size(): errors.append(['track_link',evt])
 physical=[]
 for p in mc:
  v=p.getMomentum()
  if not all(math.isfinite(x) for x in (v.x,v.y,v.z,p.getMass())): errors.append(['truth_nonfinite',evt])
  for q in list(p.getParents())+list(p.getDaughters()):
   if not 0<=q.getObjectID().index<mc.size(): errors.append(['truth_link',evt])
  if abs(p.getPDG())==15:
   truthstats['tau_records']+=1
   if any(abs(q.getPDG())==16 for q in p.getDaughters()):
    physical.append(p.getPDG()); truthstats['physical_tau_decays']+=1
    v=p.getVertex(); e=p.getEndpoint(); truthstats['physical_tau_nonzero_flight']+=sum((getattr(e,a)-getattr(v,a))**2 for a in ('x','y','z'))>1e-12
 for l in links:
  if not (0<=l.getFrom().getObjectID().index<ps.size() and 0<=l.getTo().getObjectID().index<mc.size()): errors.append(['reco_truth_link',evt])
 # Only compare by slot where count and PDGs prove the original LU record order.
 if 1<=evt<=len(fad):
  rows=fad[evt-1]
  if len(rows)==mc.size() and all(r[1]==p.getPDG() for r,p in zip(rows,mc)):
   comparable+=1
   for r,p in zip(rows,mc):
    v=p.getMomentum(); max_pdiff=max(max_pdiff,*(abs(r[5+j]-getattr(v,a)) for j,a in enumerate(('x','y','z'))))
   # Relative vertices remove the independently smeared simulation origin.
   if mc.size():
    origin=mc[0].getVertex(); r0=rows[0]
    for r,p in zip(rows,mc):
     v=p.getVertex(); max_vdiff=max(max_vdiff,*(abs((r[10+j]-r0[10+j])-(getattr(v,a)-getattr(origin,a))) for j,a in enumerate(('x','y','z'))))
  else: truthstats['fadgen_order_or_size_diff_events']+=1
 else: errors.append(['event_outside_generator_range',run,evt])
 eventrows.append(dict(run=run,event=evt,tracks=tr.size(),charged=charged,neutral=ps.size()-charged,truth_particles=mc.size(),physical_tau_decays=len(physical),reco_truth_links=links.size()))
meta=[]
for fr in reader.get('metadata'):
 if int(fr.get_parameter('mc_weight_schema'))!=1: raise ValueError('Missing weight schema')
 ledger_ids=list(fr.get_parameter('mc_generated_event_ids'))
 ledger_weights=list(fr.get_parameter('mc_generated_weights'))
 if ledger_ids!=sorted(weight_rows) or ledger_weights!=[float(weight_rows[k]['weight']) for k in ledger_ids]: raise ValueError('Embedded generator ledger mismatch')
 if sorted(fr.get_parameter('mc_converted_event_ids'))!=sorted(e for r,e in ids): raise ValueError('Embedded converted IDs mismatch')
 for key,value in [('mc_converted_sumw',math.fsum(observed_weights)),('mc_converted_sumw2',math.fsum(w*w for w in observed_weights))]:
  if not math.isclose(float(fr.get_parameter(key)),value,rel_tol=1e-10,abs_tol=1e-10): raise ValueError('Embedded weight sum mismatch')
 meta.append({key:fr.get_parameter(key) for key in ('dst_pa_modules_present','dst_pilot_blocklets_present','skelana_IFLCUT','skelana_IFLSTR')})
if len(meta)!=1: raise ValueError('Expected one production metadata frame')
summary=dict(events=len(ids),unique_identities=len(set(ids)),run_counts=dict(collections.Counter(r for r,e in ids)),event_min=min((e for r,e in ids),default=None),event_max=max((e for r,e in ids),default=None),ecms_counts=dict(ecs),processing_tags=dict(tags),collections=sizes,metadata=meta,truth_stats=dict(truthstats),fadgen_slot_comparable_events=comparable,max_momentum_difference_gev=max_pdiff,max_relative_vertex_difference_mm=max_vdiff,errors=errors)
summary['weight_audit']=dict(passed=True,count=len(observed_weights),sumw=math.fsum(observed_weights),sumw2=math.fsum(w*w for w in observed_weights),negative_count=sum(w<0 for w in observed_weights),min=min(observed_weights,default=None),max=max(observed_weights,default=None))
summary['reco_diagnostics']=reco_diagnostics
if len(ids)!=len(set(ids)): errors.append(['duplicate_identities'])
(outdir/'root-audit.json').write_text(json.dumps(summary,indent=2,default=str)+'\n')
with (outdir/'converted-events.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(eventrows[0]) if eventrows else ['run','event']); w.writeheader(); w.writerows(eventrows)
print(json.dumps({k:v for k,v in summary.items() if k!='collections'},indent=2,default=str))
sys.exit(bool(errors) or not ids)

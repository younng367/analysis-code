from pathlib import Path
import csv,json,sys
root=Path(sys.argv[1]);ranks=list(csv.DictReader((root/'data/network/candidate_rankings.csv').open(encoding='utf-8-sig')))
assert len(ranks)==48
selected={}
for alg in ['DMNC','MCC','MNC']:
 with (root/'data/network'/f'{alg.lower()}.csv').open() as f:next(f);rows=list(csv.DictReader(f))
 selected[alg]=[r['Name'] for r in rows[:3]]
 subset=[r for r in ranks if r['Algorithm']==alg];assert len(subset)==16
 for position,row in enumerate(rows,1):
  item=next(r for r in subset if r['Network_label']==row['Name'])
  assert item['Score']==row['Score'] and int(item['Tied_rank'])==int(row['Rank']) and int(item['Export_position'])==position
 assert sum(r['In_archived_score_export']=='False' and r['Score']=='NA' for r in subset)==6
intersection=sorted(set.intersection(*(set(s) for s in selected.values())))
assert intersection==['Coro1a','Lcp1']
print(json.dumps({'fixed_top3':selected,'intersection':intersection,'verified_exported_scores':30,'table_records':48},indent=2))

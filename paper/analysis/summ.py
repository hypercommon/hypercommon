import json,sys
for n in sys.argv[1:]:
    d=json.load(open(f'{n}.json'));e=d['edges'];w=d['wedges']
    print('=====',n,e['class_counts'],'AUC in/out',round(e['auc_inner_vs_outer'],3),'pooled',round(e['auc_inner_vs_pooled'],3),'base',round(e['base_rate_inner'],3))
    print({k:(v and round(v['median'],3)) for k,v in e['jaccard'].items()})
    for s in e['degree_stratified']: print('  deg',s['dmax_lo'],s['dmax_hi'],s['n_inner'],s['n_outer'],round(s['median_inner'],3),round(s['median_outer'],3),'auc',round(s['auc'],3))
    print('no-tri',e['frac_edges_no_triangle'])
    print('wedge base',round(w['base_rate_pure'],3),'n',w['evaluated_labelled'],'auc all',round(w['auc_min_all'],3),'auc edges',round(w['auc_min_edges'],3),'closed pure/mixed',w['frac_closed_pure'],w['frac_closed_mixed'])
    for c,p in zip(w['curve'][:12],e['pair_rule'][:12]): print('  t',c['t'],'edge-prec',round(p['precision'],3),'edge-rec',round(p['recall'],3),'| wedge pair-prec',round(c['pair_rule_precision'],3),'triple-prec',round(c['triple_rule_precision'],3),'triple-rec',round(c['triple_rule_recall'],3))

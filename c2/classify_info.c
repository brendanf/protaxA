#include <time.h>
#include "defs.h"

int compute_cnode_probs(TaxonomyNode *node, int nid, double prevprob, Model *m, double **scs, double pth, double *pdistances) {
  int i,j,cid,k;
  double dist, mindist1, mindist2, maxz,ezsum, *beta, *sc;

  beta = m->params[node[nid].level];
  sc = scs[node[nid].level];
  maxz = 0.0;
  for (i=0; i<node[nid].num_cnodes; i++) {
    cid = node[nid].cnode_index[i];
    mindist1 = 1.0;
    mindist2 = 1.0;
    node[cid].ind1 = -1;
    node[cid].ind2 = -1;
    for (j=0; j<node[cid].num_rseqs; j++) {
      k = node[cid].rseq_index[j];
      /* dist = pdist(seq, rseq->seq[k], rseq->alen); */
      dist = pdistances[k];
      if (dist < mindist1) {
	mindist2 = mindist1;
	mindist1 = dist;
	node[cid].dist2 = mindist2;
	node[cid].dist1 = mindist1;
	node[cid].ind2 = node[cid].ind1;
	node[cid].ind1 = k;
      }
      else if (dist < mindist2) {
	mindist2 = dist;
	node[cid].dist2 = mindist2;
	node[cid].ind2 = k;
      }
    }

    /* printf("  %s %f %f\n",node[cid].name,mindist,avedist); */

    /* use prob temporarily to store z */
    if (node[cid].isunk) {
      node[cid].prob = 0.0;
      node[cid].no_rseqs = 1;
    }
    else if (node[cid].num_rseqs) {
      if (node[cid].num_rseqs==1) mindist2=mindist1;
      node[cid].prob = beta[1] + beta[2]*(mindist1-sc[0])/sc[1] + beta[3]*(mindist2 - mindist1 - sc[2])/sc[3];
      node[cid].no_rseqs = 0;
    }
    else {
      node[cid].prob = beta[0];
      node[cid].no_rseqs = 1;
    }

    if (node[cid].prob > maxz)
      maxz = node[cid].prob;
  }

  ezsum = 1e-100;
  for (i=0; i<node[nid].num_cnodes; i++) {
    cid = node[nid].cnode_index[i];
    /* node[cid].prob = node[cid].prior * exp(node[cid].prob - maxz); */
    node[cid].prob = exp(node[cid].prob - maxz);
    ezsum += node[cid].prob;
  }

  for (i=0; i<node[nid].num_cnodes; i++) {
    cid = node[nid].cnode_index[i];
    node[cid].prob /= ezsum;
  }

  node[nid].sumcprob_no_rseqs = 0.0;
  for (i=0; i<node[nid].num_cnodes; i++) {
    cid = node[nid].cnode_index[i];
    node[cid].prob *= prevprob;
    if (node[cid].no_rseqs)
      node[nid].sumcprob_no_rseqs += node[cid].prob;
    if ((node[cid].no_rseqs == 0) && (node[cid].prob >= pth)) {
      printf("%d %s %f dist1: %d %.3f dist2: %d %.3f \n",node[cid].level,node[cid].name, node[cid].prob, node[cid].ind1, node[cid].dist1,node[cid].ind2, node[cid].dist2);
      if (node[cid].num_cnodes)
	compute_cnode_probs(node, cid, node[cid].prob, m, scs, pth, pdistances);
    }
  }
  if (node[nid].sumcprob_no_rseqs >= pth) {
    if (node[nid].level == 0)
      printf("%d %s %f\n",node[nid].level+1, UNKNAME, node[nid].sumcprob_no_rseqs);
    else
      printf("%d %s,%s %f\n",node[nid].level+1,node[nid].name, UNKNAME, node[nid].sumcprob_no_rseqs);
  }

  return(0);
}


int main (int argc, char **argv) {
  InputOptions iopt;
  char *rfile, *ifile;
  int i,j, num_tnodes, num_sclevels;
  SequenceSetB *rseq,*iseq;
  TaxonomyNode *taxonomy;
  Model *model;
  double pth, **scs;
  double *pdistances;

  iopt = get_input_options(argc, argv);

  if (argc - optind != 7) {
    fprintf(stderr,"classify sequences and print out information regarding distances to refseqs\n");
    fprintf(stderr,"usage: classify_info [-l len] [-r n_rseq] [-i n_iseq] taxonomy rseqFASTA taxonomy2rseq modelparameters scalingfile probability_threshold inputFASTA\n");
    exit(0);
  }

  taxonomy = read_taxonomy(argv[optind++], &num_tnodes);
  rfile = argv[optind++];
  add_rseq2taxonomy(argv[optind++], taxonomy);
  model = read_model(argv[optind++]);
  scs=read_level_scalings(argv[optind++], &num_sclevels);

  if (model->num_levels != num_sclevels) {
    fprintf(stderr,"ERROR: %d model levels but %d scaling levels, files '%s' and '%s'.\n", model->num_levels, num_sclevels, argv[4], argv[5]);
    exit(0);
  }

  pth = atof(argv[optind++]);
  ifile = argv[optind++];

  read_sequence_setsB(iopt, rfile, ifile, &rseq, &iseq);

  /*
    print_taxonomy(taxonomy, num_tnodes);
    print_model(model);
  */

  if ((pdistances = (double *) malloc(rseq->num_seqs * sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR: cannot malloc %d doubles for pdistances.\n",rseq->num_seqs);
    perror(""); exit(-1);
  }

  for (i=0; i<iseq->num_seqs; i++) {
    compute_distancesB(rseq, iseq->b[i], iseq->m[i], pdistances);
    printf("%s\n",iseq->id[i]);
    compute_cnode_probs(taxonomy, 0, 1.0, model, scs, pth, pdistances);
    printf("\n");
  }

  return(0);
}

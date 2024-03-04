#include "defs.h"

int *malloc_int_array(int n) {
  int *a;
  if ((a = (int *) malloc (n * sizeof(int))) == NULL) {
    fprintf(stderr,"ERROR: cannot malloc %d ints.\n",n);
    perror(""); exit(-1);
  }
  return(a);
}

int *read_rseq_indices(char *filename, int *num_indices) {
  FILE *fp;
  char line[MAXLINE], *token;
  int *rindex;
  int i, j, n;
  
  if ((fp = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"ERROR: cannot read input rseq indices '%s'.\n",filename);
    perror(""); exit(-1);
  }

  n=0;
  while (fgets(line,MAXLINE,fp) != NULL) {
    n++;
  }

  rindex = malloc_int_array(n);
  
  rewind(fp);

  for (j=0; j<n; j++) {
    if (fgets(line,MAXLINE,fp) == NULL) {
      fprintf(stderr,"ERROR: cannot read index %d/%d from '%s'.\n",j+1,n,filename);
      perror(""); exit(-1);
    }
      
    if ((token = strtok(line," \t")) == NULL) {
      fprintf(stderr,"ERROR: cannot read %d/%d token from file '%s'.\n",j+1,n,filename);
      perror(""); exit(-1);
    }
    rindex[j] = atoi(token);
  }
  fclose(fp);
  
  *num_indices = n;
  return (rindex);
}

int compute_cnode_probs(TaxonomyNode *node, int nid, double prevprob, Model *m, double **scs, double pth, double *pdistances, int self_index) {
  int i,j,cid,k;
  double dist,mindist1, mindist2, maxz,ezsum, *beta, *sc;
  int num_rseqs;

  beta = m->params[node[nid].level];
  sc = scs[node[nid].level];
  maxz = 0.0;
  for (i=0; i<node[nid].num_cnodes; i++) {
    cid = node[nid].cnode_index[i];
    mindist1 = 1.0;
    mindist2 = 1.0;
    num_rseqs = node[cid].num_rseqs;
    for (j=0; j<node[cid].num_rseqs; j++) {
      k = node[cid].rseq_index[j];
      if (k == self_index) {
	num_rseqs--;
      }
      else {
	/* dist = pdist(seq, rseq->seq[k], rseq->alen); */
	dist = pdistances[k];
	if (dist < mindist1) {
	  mindist2 = mindist1;
	  mindist1 = dist;
	}
	else if (dist < mindist2) {
	  mindist2 = dist;
	}
      }
    }
    
    /* printf("  %s %f %f\n",node[cid].name,mindist,avedist); */
    
    /* use prob temporarily to store z */
    if (node[cid].isunk) {     
      node[cid].prob = 0.0;
      node[cid].no_rseqs = 1;
    }
    else if (num_rseqs) {
      if (num_rseqs==1) mindist2=mindist1;
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
    /* node[cid].prob = node[cid].prior * exp(node[cid].prob - maxz);*/
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
      printf(" %s %f",node[cid].name, node[cid].prob);
      if (node[cid].num_cnodes)
	compute_cnode_probs(node, cid, node[cid].prob, m, scs, pth, pdistances, self_index);
    }
  }
  if (node[nid].sumcprob_no_rseqs >= pth) {
    if (node[nid].level == 0)
      printf(" %s %f",UNKNAME, node[nid].sumcprob_no_rseqs);
    else
      printf(" %s,%s %f",node[nid].name, UNKNAME, node[nid].sumcprob_no_rseqs);
  }

  return(0);
}


int main (int argc, char **argv) {
  int i,j, num_tnodes, num_sclevels, rlen, n_rseq;
  SequenceSetB *rseq;
  TaxonomyNode *taxonomy;
  Model *model;
  double pth, **scs;
  double *pdistances;
  int n_input_index, *input_index;

  if (argc < 7) {
    fprintf(stderr, "classify reference sequences (whose indices in index_file) by not utilizing self-similarity\n");
    fprintf(stderr,"usage: classify_rseq taxonomy rseqFASTA taxonomy2rseq modelparameters scalingfile probability_threshold input_rindex_file\n");
    exit(0);	    
  }

  taxonomy = read_taxonomy(argv[1], &num_tnodes);
  scan_aligned_sequences(argv[2], &rlen, &n_rseq);
  rseq = read_aligned_sequencesB(argv[2], rlen, n_rseq);
  add_rseq2taxonomy(argv[3], taxonomy);
  model = read_model(argv[4]);
  scs=read_level_scalings(argv[5], &num_sclevels);

  if (model->num_levels != num_sclevels) {
    fprintf(stderr,"ERROR: %d model levels but %d scaling levels, files '%s' and '%s'.\n", model->num_levels, num_sclevels, argv[4], argv[5]);
    exit(0);
  }
  
  pth = atof(argv[6]);  
  input_index = read_rseq_indices(argv[7], &n_input_index);
  
  if ((pdistances = (double *) malloc(rseq->num_seqs * sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR: cannot maloc %d doubles for pdistances.\n",rseq->num_seqs);
    perror(""); exit(-1);
  }

  for (i=0; i<n_input_index; i++) {
    j = input_index[i];
    compute_distancesB(rseq, rseq->b[j], rseq->m[j], pdistances);
    printf("%s",rseq->id[j]);
    compute_cnode_probs(taxonomy, 0, 1.0, model, scs, pth, pdistances, j);
    printf("\n");
  }
  
  return(0);
}

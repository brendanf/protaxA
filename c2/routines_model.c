#include "defs.h"
#include <math.h>

Model *read_model(char *filename) {
  FILE *fp;
  char line[MAXLINE], *token;
  int i,j;
  Model *m;
  char *thisfunction = "read_model";

  if ((fp = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot open '%s' for reading.\n",thisfunction,filename);
    perror(""); exit(-1);
  }

  if ((m = (Model *) malloc (sizeof(Model))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc Model.\n", thisfunction);
    perror(""); exit(-1);
  }

  m->num_levels=0;
  m->dim = 0;

  while (fgets(line,MAXLINE,fp))
    m->num_levels++;

  if (strtok(line," \t\n") == NULL) {
    fprintf(stderr,"ERROR (%s): empty line in '%s'.\n",thisfunction, filename);
    exit(-1);
  }
  m->dim=1;
  while (strtok(NULL," \t\n"))
    m->dim++;

  if ((m->params = (double **) malloc (m->num_levels * sizeof(double *))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d*%d doubles, file '%s'.\n",thisfunction,m->num_levels,m->dim,filename);
    perror(""); exit(-1);
  }
  for (i=0; i<m->num_levels; i++) {
    if ((m->params[i] = (double *) malloc (m->dim * sizeof(double))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d*%d doubles, file '%s'.\n",thisfunction,m->num_levels,m->dim,filename);
      perror(""); exit(-1);
    }
  }

  rewind(fp);

  for (i=0; i<m->num_levels; i++) {
    if (fgets(line,MAXLINE,fp) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot re-read level %d/%d parameters, file '%s'.\n",thisfunction,i+1,m->num_levels,filename);
      perror(""); exit(-1);
    }
    if ((token = strtok(line," \t\n")) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot re-read level %d/%d 1st parameter, file '%s'.\n",thisfunction,i+1,m->num_levels,filename);
      perror(""); exit(-1);
    }
    j=0;
    m->params[i][j] = atof(token);
    for (j=1; j<m->dim; j++) {
      if ((token = strtok(NULL," \t\n")) == NULL) {
	fprintf(stderr,"ERROR (%s): cannot re-read level %d/%d parameter %d/%d, file '%s'.\n",thisfunction,i+1,m->num_levels,j+1,m->dim,filename);
	perror(""); exit(-1);
      }
      m->params[i][j] = atof(token);
    }
  }
  fclose(fp);

  return(m);
}

double **read_level_scalings(char *filename, int *num_levels) {
  FILE *fp;
  char line[MAXLINE], *token;
  double **sc;
  int i, j, nlev;

  if ((fp = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"ERROR: cannot read scalingfile '%s'.\n",filename);
    perror(""); exit(-1);
  }

  nlev=0;
  while (fgets(line,MAXLINE,fp) != NULL) {
    nlev++;
  }

  if ((sc = (double **) malloc (nlev*sizeof(double *))) == NULL) {
    fprintf(stderr,"ERROR: cannot malloc %d double ptrs.\n", nlev);
    perror(""); exit(-1);
  }
  for (j=0; j<nlev; j++) {
    if ((sc[j] = (double *) malloc (4*sizeof(double))) == NULL) {
      fprintf(stderr,"ERROR: cannot malloc %d*4 doubles.\n",nlev);
      perror(""); exit(-1);
    }
  }

  rewind(fp);

  for (j=0; j<nlev; j++) {
    if (fgets(line,MAXLINE,fp) == NULL) {
      fprintf(stderr,"ERROR: cannot read scaling factors %d/%d from '%s'.\n",j+1,nlev,filename);
      perror(""); exit(-1);
    }

    if ((token = strtok(line," \t")) == NULL) {
      fprintf(stderr,"ERROR: cannot read level %d 1st token from file '%s'.\n",nlev,filename);
      perror(""); exit(-1);
    }
    if ((token = strtok(NULL," \t")) == NULL) {
      fprintf(stderr,"ERROR: cannot read level %d 2nd token from file '%s'.\n",nlev,filename);
      perror(""); exit(-1);
    }
    sc[j][0] = atof(token);
    for (i=1; i<4; i++) {
      if ((token = strtok(NULL," \t")) == NULL) {
	fprintf(stderr,"ERROR: cannot read level %d token %d from file '%s'.\n",nlev,2*i+1,filename);
	perror(""); exit(-1);
      }
      if ((token = strtok(NULL," \t")) == NULL) {
	fprintf(stderr,"ERROR: cannot read level %d token %d from file '%s'.\n",nlev,2*i+2,filename);
	perror(""); exit(-1);
      }
      sc[j][i] = atof(token);
    }
  }
  fclose(fp);

  *num_levels = nlev;
  return (sc);
}

int compute_cnode_probs_best2(const char *qid, TaxonomyNode *node, int nid, double prevprob, const Model *m, const double **scs, double pth, double rth, const double *pdistances) {
  int i,j,cid,k;
  double dist,mindist1, mindist2, maxz,ezsum, *beta;
  const double *sc;

  /* model utilizes distances to best and 2nd best refseq neighbors in each node (mindist1, mindist2) */

  beta = m->params[node[nid].level];
  sc = scs[node[nid].level];
  maxz = 0.0;
  for (i=0; i<node[nid].num_cnodes; i++) {
    cid = node[nid].cnode_index[i];
    mindist1 = 1.0;
    mindist2 = 1.0;
    for (j=0; j<node[cid].num_rseqs; j++) {
      k = node[cid].rseq_index[j];
      dist = pdistances[k];
      if (dist < mindist1) {
	mindist2 = mindist1;
	mindist1 = dist;
      }
      else if (dist < mindist2) {
	mindist2 = dist;
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
    if (node[cid].no_rseqs == 0) {
      if (node[cid].prob >= rth) {
        printf("%s %d %s %f\n", qid, node[cid].level, node[cid].name, node[cid].prob);
      }
      if (node[cid].num_cnodes && (node[cid].prob >= pth))
        compute_cnode_probs_best2(qid, node, cid, node[cid].prob, m, scs, pth, rth, pdistances);
    }
  }
  if (node[nid].sumcprob_no_rseqs >= rth) {
    if (node[nid].level == 0)
      printf("%s %d %s %f\n",qid, node[nid].level + 1, UNKNAME, node[nid].sumcprob_no_rseqs);
    else
      printf("%s %d %s,%s %f\n",qid, node[nid].level + 1, node[nid].name, UNKNAME, node[nid].sumcprob_no_rseqs);
  }

  return(0);
}

int print_model(Model *m) {
  int i,j;
  for (i=0; i<m->num_levels; i++) {
    printf("model level %d:",i+1);
    for (j=0; j<m->dim; j++)
      printf(" %f",m->params[i][j]);
    printf("\n");
  }
  return(0);
}


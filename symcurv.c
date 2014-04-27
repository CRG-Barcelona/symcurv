#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* symcurv - adaptation of algorithm and code from perl script */
/* originally by Christoforos Nikolau. The program is designed */
/* to find symmetry in curvature of DNA based on thermodynamic */
/* parameters of DNA base pairing/stacking.                    */

#include <jkweb/common.h>
#include <jkweb/linefile.h>
#include <jkweb/hash.h>
#include <jkweb/options.h>
#include <jkweb/sqlNum.h>
#include <jkweb/dnaLoad.h>
#include <jkweb/basicBed.h>
#include <jkweb/rangeTree.h>
#include <beato/bigs.h>

/* local headers: */

#include "symcurv_constants.h"

void usage()
/* Explain usage and exit. */
{
errAbort(
  "symcurv - run SymCurv symmetry-of-curvature algorithm on DNA (twoBit/FASTA)\n"
  "usage:\n"
  "   symcurv sequence_file:range output.wig\n"
  "options:\n"
  "   -regions=bed            specify specific regions to run algorithm (BED3)\n"
/*   "   -track-line=\"track...\"  add a custom-track line to the top of the output\n" */
  "   -decimals=n             specify decimal precision of output\n"
  "   -verbose                log some runtime information to stdout\n"
  "   -wigtype=bg,fix,var     type of wiggle output (default is fixedStep for\n"
  "                           curvature-only output, bedGraph for SymCurv output\n"
  "   -with-zeros             in the case that -wig=type=bg (bedGraph), then choose\n"
  "                           to output all the zeros for SymCurv output. Normally,\n"
  "                           they are skipped\n"
  "   -curvature              output curvature values as opposed to symmetry-of-\n"
  "                           curvature values.\n"
  "   -averaged               output averaged symmetry-of-curvature values at each\n"
  "                           base.\n"
  "   -nuc-calls=file.bg      output non-overlapping nucleosome calls in bedGraph format\n"
  "   -activated              use \"activated\" set of DNA parameters\n"
  "   -N-to-A                 instead of the default behavior of skipping N bases,\n"
  "                           convert N bases to A.\n"
  "   -min-length=n           run SymCurv on sequences longer than n bp (default 1000)\n"
    );
}

static struct optionSpec options[] = {
   {"regions", OPTION_STRING},
   {"andy", OPTION_BOOLEAN},
   {"wigtype", OPTION_STRING},
   {"with-zeros", OPTION_BOOLEAN},
   {"track-line", OPTION_STRING},
   {"decimals", OPTION_INT},
   {"verbose", OPTION_BOOLEAN},
   {"averaged", OPTION_BOOLEAN},
   {"N-to-A", OPTION_BOOLEAN},
   {"activated", OPTION_BOOLEAN},
   {"curvature", OPTION_BOOLEAN},
   {"min-length", OPTION_INT},
   {"nuc-calls", OPTION_STRING},
   {NULL, 0},
};

void delete_short_subsections(struct bed **pBedList)
/* delete all the subsections shorter than 500 or -min-length bases */
{
    struct bed *section = NULL;
    struct bed *newList = NULL;
    int min_len = optionInt("min-length", 1000);
    if (min_len < 2 * SC_SYMCURVE_WIN + 3)
	errAbort("error: -min-length must be >= %d", 2 * SC_SYMCURVE_WIN + 3);
    while ((section = slPopHead(pBedList)) != NULL)
	if (section->chromStart - section->chromEnd >= min_len)
	    slAddHead(&newList, section);
	else
	    bedFree(&section);
    slReverse(&newList);
    *pBedList = newList;
}

void transform_seq(struct dnaSeq *seq)
/* convert sequence to all upper case, and if necessary, convert N->A */
{
    int i;
    touppers(seq->dna);
    if (optionExists("N-to-A"))
	for (i = 0; i < seq->size; i++)
	    if (seq->dna[i] == 'N')
		seq->dna[i] = 'A';
}

int *get_num_seq(struct dnaSeq *seq, struct bed *section)
/* translate DNA sequence into indeces into roll/tilt/whatever arrays */
/* that expect A=0, T=1, G=2, and C=3 */
{
    int i;
    int *ret;
    int size = section->chromEnd - section->chromStart;
    AllocArray(ret, section->chromEnd - section->chromStart);
    for (i = 0; i < size; i++)
    {
	switch (seq->dna[i + section->chromStart])
	{
	case 'A':
	    ret[i] = 0;
	    break;
	case 'T':
	    ret[i] = 1;
	    break;
	case 'G':
	    ret[i] = 2;
	    break;
	case 'C':
	    ret[i] = 3;
	    break;
	default:
	    errAbort("Non-DNA base found in input sequence");
	}
    }
    return ret;
}

void do_curv(struct perBaseWig *curv, struct dnaSeq *seq, struct bed *section)
/* Do the actual curvature calculations */
{
    int i = 0;
    int size = section->chromEnd - section->chromStart;
    int *ns = get_num_seq(seq, section);
    double twist_sum = 0;
    double *xcoord = NULL;
    double *ycoord = NULL;
    double *xave = NULL;
    double *yave = NULL;
    boolean use_dnase_roll = optionExists("activated");
    AllocArray(xcoord, size);
    AllocArray(ycoord, size);
    AllocArray(xave, size);
    AllocArray(yave, size);
    /* first pass... */
    for (i = 0; i < size - 2; i++)
    {
	double roll = 0;
	double dx, dy;
	twist_sum += SC_TWIST[ns[i]][ns[i+1]][ns[i+2]];
	if (use_dnase_roll)
	    roll = SC_ROLL_DNASE[ns[i]][ns[i+1]][ns[i+2]];
	else
	    roll = SC_ROLL_NUC[ns[i]][ns[i+1]][ns[i+2]];
	dx = roll * sin(twist_sum);
	dy = roll * cos(twist_sum);
	xcoord[i+1] = xcoord[i] + dx;
	ycoord[i+1] = ycoord[i] + dy;
    }
    /* second pass */
    for (i = SC_CURVE_STEPONE; i < size - SC_CURVE_STEPONE; i++)
    {
	double rxsum = 0;
	double rysum = 0;
	int j;
	for (j = -1 * SC_CURVE_STEPTWO; j <= SC_CURVE_STEPTWO; j++)
	{
	    rxsum += xcoord[i+j];
	    rysum += ycoord[i+j];
	}
	rxsum += xcoord[i+(SC_CURVE_STEPONE-1)]/(SC_CURVE_STEPTWO/2) + 
	    xcoord[i-(SC_CURVE_STEPONE-1)]/(SC_CURVE_STEPTWO/2);
	rysum += ycoord[i+(SC_CURVE_STEPONE-1)]/(SC_CURVE_STEPTWO/2) + 
	    ycoord[i-(SC_CURVE_STEPONE-1)]/(SC_CURVE_STEPTWO/2);
	xave[i] = rxsum/((SC_CURVE_STEPONE-1)*2);
	yave[i] = rysum/((SC_CURVE_STEPONE-1)*2);
    }
    /* third pass: fill in the wig */
    for (i = SC_CURVE_STEP + SC_CURVE_STEPONE; i < size - SC_CURVE_STEP - SC_CURVE_STEPONE; i++)
    {
	double x_square_part = pow(xave[i+SC_CURVE_STEP] - xave[i-SC_CURVE_STEP], 2); 
	double y_square_part = pow(yave[i+SC_CURVE_STEP] - yave[i-SC_CURVE_STEP], 2); 
	curv->data[section->chromStart + i] = sqrt(x_square_part + y_square_part) * SC_CURVE_SCALE;
    }
    freeMem(xcoord);
    freeMem(ycoord);
    freeMem(xave);
    freeMem(yave);
    freeMem(ns);
}

double do_symcurv(struct perBaseWig *curv, struct perBaseWig *symcurv, struct bed *section)
/* the final calculation: the symmetry-of-curvature. return the average across this window */
{
    int win = SC_SYMCURVE_WIN;
    int step = SC_SYMCURVE_STEP;
    int dyad;
    int size = section->chromEnd - section->chromStart;
    double *curv_nuc = curv->data + section->chromStart;
    double *symcurv_data = symcurv->data + section->chromStart;
    double average = 0;
    /* calculate across region */
    for (dyad = win; dyad < size - win; dyad += step)
    {
	double weight = 0;
	double sum = 0; 
	int j, k;
	double inv_weight = (curv_nuc[dyad-1] - curv_nuc[dyad]) + (curv_nuc[dyad+1] - curv_nuc[dyad]);
	for (j = dyad, k = dyad; (j < dyad + win/2 + 1) && (k > (dyad - win/2 - 1)); j += step, k -= step)
	    sum += fabs(curv_nuc[j]-curv_nuc[k]);
	/* Right, so this is the interesting bit that zeroes-out the score if specific restrictions aren't met: */
	/*  (1)  The curvature at the dyad must be less than the base before and the base after. */
	/*  (2)  That combined difference is greater than 0.01. */
	if ((curv_nuc[dyad] < curv_nuc[dyad-1]) && (curv_nuc[dyad] < curv_nuc[dyad+1]) && (inv_weight >= 0.01))
	    weight = 1/inv_weight;
	else
	    weight = 0;
	if (sum != 0)
	    symcurv_data[dyad] = 1/sum * weight;
	else 
	    symcurv_data[dyad] = 100;
	average += symcurv_data[dyad];
    }
    if (size - 2*win > 0)
	average /= size - 2*win;
    return average;
}

struct spear_rank
{
    int coord;
    double x_i;
    double y_i;
    double val;
};

int spear_rank_cmp(const void *va, const void *vb)
{
    struct spear_rank *a = (struct spear_rank *)va;
    struct spear_rank *b = (struct spear_rank *)vb; 
    if (a->val < b->val)
	return -1;
    if (b->val < a->val)
	return 1;
    return 0;
}

void dump_corr_ranks(struct spear_rank corr_ranks[7])
{
    int i;
    for (i = 0; i < 7; i++)
    {
	uglyf("cr[%d].val = %f\tcr[%i].coord = %d\tcr[%d].x_i = %0.2f\tcr[%d].y_i = %0.2f\n", i, corr_ranks[i].val, i, corr_ranks[i].coord, i, corr_ranks[i].x_i, i, corr_ranks[i].y_i);
    }
}

double compute_spearman(struct spear_rank corr_ranks[7])
{
    int i=0, j, k;
    int sum;
    double mean_x_i = 0;
    double mean_y_i = 0;
    double numer = 0;
    double denom = 0;
    double denom_sumx = 0;
    double denom_sumy = 0;
    qsort(corr_ranks, 7, sizeof(struct spear_rank), spear_rank_cmp);
    /* fill in x_i values and adjust for tied-ranks */
    while (i < 7)
    {
	j = i+1;
	sum = i+1;
	while ((j < 7) && (corr_ranks[i].val == corr_ranks[j].val))
	{
	    j++;
	    sum += j;
	}
	if ((j-i) > 1)
	{
	    for (k = 0; k < (j-i); k++)
		corr_ranks[i+k].x_i = (double)sum/(j-i);
	}
	else /* no ties */
	    corr_ranks[i].x_i = (double)i+1;
	i = j;
    }
    /* calculate x_i/y_i means */
    for (i = 0; i < 7; i++)
    {
	mean_x_i += corr_ranks[i].x_i;
	mean_y_i += corr_ranks[i].y_i;
    }
    mean_x_i /= 7;
    mean_y_i /= 7;
    /* replace x_i and y_i values with differences to means */
    for (i = 0; i < 7; i++)
    {
	corr_ranks[i].x_i -= mean_x_i;
	corr_ranks[i].y_i -= mean_y_i;	
    }
    /* one more time go through and do the top and bottom of the final fraction */
    for (i = 0; i < 7; i++)
    {
	numer += corr_ranks[i].x_i * corr_ranks[i].y_i;
	denom_sumx += pow(corr_ranks[i].x_i,2);
	denom_sumy += pow(corr_ranks[i].y_i,2);	
    }
    denom = sqrt(denom_sumx * denom_sumy);
    return numer / denom;
}

void init_corr_ranks(struct spear_rank corr_ranks[7])
{
    corr_ranks[0].coord = -100;
    corr_ranks[0].x_i = 0;
    corr_ranks[0].y_i = 4.5;
    corr_ranks[0].val = 0;
    corr_ranks[1].coord = -80;
    corr_ranks[1].x_i = 0;
    corr_ranks[1].y_i = 6.5;
    corr_ranks[1].val = 0;
    corr_ranks[2].coord = -35;
    corr_ranks[2].x_i = 0;
    corr_ranks[2].y_i = 1.5;
    corr_ranks[2].val = 0;
    corr_ranks[3].coord = 0;
    corr_ranks[3].x_i = 0;
    corr_ranks[3].y_i = 3;
    corr_ranks[3].val = 0;
    corr_ranks[4].coord = 35;
    corr_ranks[4].x_i = 0;
    corr_ranks[4].y_i = 1.5;
    corr_ranks[4].val = 0;
    corr_ranks[5].coord = 80;
    corr_ranks[5].x_i = 0;
    corr_ranks[5].y_i = 6.5;
    corr_ranks[5].val = 0;
    corr_ranks[6].coord = 100;
    corr_ranks[6].x_i = 0;
    corr_ranks[6].y_i = 4.5;
    corr_ranks[6].val = 0;
}

void do_andy_symcurv(struct perBaseWig *curv, struct perBaseWig *symcurv)
{
    const double na = NANUM;
    /* int max = (sym_ranges->max > corr_ranges->max) ? sym_ranges->max : corr_ranges->max; */
    int max = 50;
    int edge_size = max + SC_CURVE_STEP + SC_CURVE_STEPONE;
    int i, j, k;
    struct spear_rank corr_ranks[7];
    if (max < symcurv->len)
    {
	for (i = 0; i < edge_size; i++)
	    symcurv->data[i] = na;
	for (i = symcurv->len - edge_size; i < symcurv->len; i++)
	    symcurv->data[i] = na;    
	
	for (i = edge_size; i < symcurv->len - edge_size; i++)
	{
	    double dist = 0;
	    double corr = 0;
            /* euclidean distance -/+ 100 bases */
	    for (k = 1; k <= 100; k++)
		dist += pow(curv->data[i-k] - curv->data[i+k],2);
	    dist = sqrt(dist);
	    dist = 1 - dist/201;
	    /* in this case, just points are used because curvature is calculated over */
	    /* a window anyway. */
	    /* for now this is fixed at bases -/+100, -/+80, -/+35, and 0. */
	    /* ... but first, init array */
	    init_corr_ranks(corr_ranks);
	    /* load up the values */
	    corr_ranks[0].val = curv->data[i-100]; 
	    corr_ranks[1].val = curv->data[i-80];
	    corr_ranks[2].val = curv->data[i-35]; 
	    corr_ranks[3].val = curv->data[i];
	    corr_ranks[4].val = curv->data[i+35]; 
	    corr_ranks[5].val = curv->data[i+80];
	    corr_ranks[6].val = curv->data[i+100]; 
	    corr = compute_spearman(corr_ranks);
	    if (corr < 0)
		corr = 0;
	    symcurv->data[i] = corr * dist;
	}
    }
    else
    {
	for (i = 0; i < symcurv->len; i++)
	    symcurv->data[i] = na;
    }
}

struct nuc_call
/* just a simple struct to hold the relevant data */
/* that needs sorting by score */
{
    struct nuc_call *next;
    int start;
    double score;
};

int nuc_call_Cmp(const void *va, const void *vb)
/* sort callback function for nucleosome calls */
{
    const struct nuc_call *a = *((struct nuc_call **)va);
    const struct nuc_call *b = *((struct nuc_call **)vb);
    double d = b->score - a->score;;
    if (d < 0)
	return -1;
    if (d > 0)
	return 1;
    return 0;
}

struct nuc_call *nonzero_pbw(struct perBaseWig *sc)
/* return all overlapping nucleosome calls (nonzero symcurv values) */
/* sorted by score in descending order */
{
    int size = sc->chromEnd - sc->chromStart;
    int i;
    struct nuc_call *list = NULL;
    for (i = SC_SYMCURVE_WIN; i < size - SC_SYMCURVE_WIN; i++)
	if (sc->data[i] > 0)
	{
	    struct nuc_call *nc;
	    AllocVar(nc);
	    nc->start = i;
	    nc->score = sc->data[i];
	    slAddHead(&list, nc);
	}
    slSort(&list, nuc_call_Cmp);
    return list;
}
	
struct bed *nucleosome_calls(struct perBaseWig *sc, int decimals)
/* greedy algorithm to determine nucleosome calls */
{
    struct nuc_call *nuc_calls = nonzero_pbw(sc);
    struct nuc_call *call;
    struct rbTree *tree = rangeTreeNew();
    struct bed *nuc_list = NULL;
    struct range *range_list = NULL;
    struct range *range;
    /* put all the nucleosome call ranges in a red-black tree */
    for (call = nuc_calls; call != NULL; call = call->next)
    {
	int start = call->start - 74;
	int end = start + 147;
	int spacer = SC_MIN_LINKER_SIZE;
	if (!rangeTreeOverlaps(tree, start - spacer, end + spacer))
	    rangeTreeAddVal(tree, start, end, &call->score, NULL);
    }
    range_list = rangeTreeList(tree);
    /* finally just convert the tree to a bedGraph */
    for (range = range_list; range != NULL; range = range->next)
    {
	struct bed *bed;
	char number[16];
	double *val = range->val;
	safef(number, sizeof(number), "%0.*f", decimals, *val);
	AllocVar(bed);
	bed->chrom = cloneString(sc->chrom);
	bed->chromStart = sc->chromStart + range->start;
	bed->chromEnd = sc->chromStart + range->end;
	bed->name = cloneString(number);
	slAddHead(&nuc_list, bed);
    }
    slReverse(&nuc_list);
    slFreeList(&nuc_calls);
    rangeTreeFree(&tree);
    return nuc_list;
}

void do_averaging(struct perBaseWig *pbw)
/* average SymCurv values over 147 bp windows */
{
    int i, j;
    double *aves;
    struct bed *section;
    AllocArray(aves, pbw->chromEnd - pbw->chromStart);
    for (section = pbw->subsections; section != NULL; section = section->next)
    {
	/* Inefficient way to average. */
	for (i = section->chromStart + 73; i < section->chromEnd - 74; i++)
	{
	    double sum = 0;
	    for (j = i - 73; j < i + 74; j++)
		sum += pbw->data[j];
	    aves[i] = sum / 147;
	}
	section->chromStart += 73;
	section->chromEnd -= 74;
    }
    freeMem(pbw->data);
    pbw->data = aves;
}

void symcurv(char *input, char *output)
/* symcurv - run SymCurv symmetry-of-curvature algorithm on DNA */
{
    int decimals = optionInt("decimals", 4);
    char *track_line = optionVal("track-line", NULL);
    boolean only_curvature = optionExists("curvature");
    boolean with_zeros = optionExists("with-zeros");
    boolean andy_symcurv = optionExists("andy");
    char *nuc_calls_file = optionVal("nuc-calls", NULL);
    boolean verbose = optionExists("verbose");
    boolean averaged = optionExists("averaged");
    struct dnaSeq *seqs = dnaLoadAll(input);
    struct dnaSeq *seq;
    enum wigOutType wot = get_wig_out_type((char *)optionVal("wigtype", "fix"));
    FILE *out = mustOpen(output, "w");
    FILE *calls_out = NULL;
    if (nuc_calls_file)
	calls_out = mustOpen(nuc_calls_file, "w");
    if (!optionExists("wigtype") && !only_curvature)
	wot = bedGraphOut;
    for (seq = seqs; seq != NULL; seq = seq->next)
    {
	struct bed *sec;
	struct perBaseWig *pbw_curv = NULL;
	struct perBaseWig *pbw_symcurv = NULL;
	int subsectionsSize = 0;
	transform_seq(seq);
	pbw_curv = alloc_perBaseWig_matchingSequence(seq, !optionExists("N-to-A"));
	delete_short_subsections(&pbw_curv->subsections);
	pbw_symcurv = perBaseWigClone(pbw_curv);
	for (sec = pbw_curv->subsections; sec != NULL; sec = sec->next)
	    subsectionsSize += sec->chromEnd - sec->chromStart;
	if (verbose)
	    printf("running symcurv on %s (%d bp) (%d bp in subsections)\n", seq->name, seq->size, subsectionsSize);
        /* running the SymCurv algorithm */
	for (sec = pbw_curv->subsections; sec != NULL; sec = sec->next)
	{
	    do_curv(pbw_curv, seq, sec);
	    if (only_curvature)
	    {
		/* remove the boundary sections of the wig where there is no data */
		sec->chromStart += SC_CURVE_STEP + SC_CURVE_STEPONE;
		sec->chromEnd -= SC_CURVE_STEP + SC_CURVE_STEPONE;
	    }
	    else
	    {
		do_symcurv(pbw_curv, pbw_symcurv, sec); 
		sec->chromStart += SC_SYMCURVE_WIN;
		sec->chromEnd -= SC_SYMCURVE_WIN;
	    }
	}
	if (averaged)
	    do_averaging(pbw_symcurv);
	if (only_curvature)
	    perBaseWigOutput(pbw_curv, out, wot, decimals, track_line, FALSE, TRUE);
	else if (andy_symcurv)
	{
	    do_andy_symcurv(pbw_curv, pbw_symcurv);
	    perBaseWigOutputNASkip(pbw_symcurv, out, wot, decimals, track_line, FALSE, TRUE);
	}
	else
	/* (symcurv) */
	{
	    if (!with_zeros && !averaged && (wot == bedGraphOut))
		perBaseWigOutputNASkip(pbw_symcurv, out, bedGraphOut, decimals, track_line, FALSE, TRUE);
	    else
		perBaseWigOutput(pbw_symcurv, out, wot, decimals, track_line, FALSE, TRUE);
	    /* do nucleosome calls */
	    if (nuc_calls_file)
	    {
		struct bed *calls = nucleosome_calls(pbw_symcurv, decimals);
		struct bed *call;
		for (call = calls; call != NULL; call = call->next)
		    bedTabOutN(call, 4, calls_out);
		bedFreeList(&calls);
	    }
	}
	perBaseWigFree(&pbw_curv);
	perBaseWigFree(&pbw_symcurv);
    }
    carefulClose(&out);
    dnaSeqFreeList(&seqs);
}

int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
if (argc != 3)
    usage();
symcurv(argv[1], argv[2]);
return 0;
}

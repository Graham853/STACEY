/*
        Copyright (C) 2015 Graham Jones, www.indriid.com

        This file is part of STACEY.

        STACEY is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        STACEY is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with STACEY.  If not, see <http://www.gnu.org/licenses/>.
*/

package stacey;

import beast.core.*;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;
import stacey.util.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 *  Created by Graham Jones on 19/08/2014.
 */

/*
This is the central class of STACEY (if any is).

It calculates the coalescent probability.
It is complicated because the calculation is quite complicated
and because it does some optimization to avoid recalculating
everything when a single gene tree is changed, and to avoid calculating
a likelihood when one or more gene trees are incompatible.
*/

/*
General notes on code organization in STACEY. 2015-07-08

Package stacey contains classes derived from BEASTObject. They all
have Inputs. These are top level classes constructed from XML.

Package stacey.util contains everything else except debugging code.
None of these classes calls methods in stacey.

Package debugtune is for debugging and tuning performance.


                                 In stacey
extend Distribution:
    BirthDeathCollapseModel, PIOMSCoalescentDistribution
extend CalculationNode:
   GtreeAndCoalFactor
extend Operator:
    CoordinatedPruneRegraft, FocusedNodeHeightScaler, NodesNudge, StaceyNodeReheight
extend BEASTObject implementing Loggable:
    BirthDeathCollapseNClustersStatistic, PopSampleStatistic
extend BEASTObject:
    InverseGammaComponent
extend InputEditor.Base
    GtreeAndCoalFactorInputEditor, PIOMSCoalescentDistributionInputEditor

The operators and this class, PIOMSCoalescentDistribution, are quite big, the rest small.

The editors aren't working.


                         In stacey.util

BitUnion is a low-level utility class providing a minimal 'Bitset' implementation.

Bindings provides various linkages between the SMC-Tree tips and gene tree tips. This is information
which is constant for the analysis. Used by operators and here.

FitsHeights is an important utility class, used for logP calculation here.
It determines whether a gene tree is compatible with the SMC-tree and
if so, assigns gene tree node heights (coalescence heights) to branches in the SMC-tree.

InverseGammaMixture implements a mixture of inverse gamma distributions.

Misc is formatting for debug output.

UnionArrays is an imporatant utility class for STACEY, used by several operators.
It adds a set of species (or minimal clusters) to every node in the SMC-tree and in all gene trees.
*/


@Description("The STACEY coalescent distribution. The central class in STACEY.")

@Citation(value = "Graham Jones (www.indriid.com). " +
        "\'STACEY: species delimitation and phylogeny estimation under the multispecies coalescent\'. " +
        "http://biorxiv.org/content/early/2015/03/22/010199",
        firstAuthorSurname = "Jones")

public class PIOMSCoalescentDistribution extends TreeDistribution {

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<List<GtreeAndCoalFactor>> geneTreesInput =
            new Input<List<GtreeAndCoalFactor>>("geneTree",
                    "All gene trees",
                    new ArrayList<GtreeAndCoalFactor>());

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<List<InverseGammaComponent>> priorComponentsInput =
            new Input<List<InverseGammaComponent>>("popPriorInvGamma",
                    "Component of mixture of inverse gamma distributions used as a prior for the per-branch populations",
                    new ArrayList<InverseGammaComponent>(),
                    Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<RealParameter> popPriorScale =
            new Input<RealParameter>("popPriorScale",
                    "Overall scale for population size",
                    Input.Validate.REQUIRED);

    // Needed for Beauti template, but not used here
    @SuppressWarnings("UnusedDeclaration")
    public Input<TaxonSet> taxonSetInput =
            new Input<TaxonSet>("taxonset",
                    "set of taxa mapping lineages to species", Input.Validate.REQUIRED);


    private Bindings bindings;
    private FitsHeights fitsHeights;
    private TreeInterface sTree;
    private int nSMCTreeNodes;
    private int nSMCTreeTips;
    private List<GtreeAndCoalFactor> gTreeCFs;
    private ArrayList<Tree> gTrees;
    private int nGTrees;

    private double [][] lnGammaRatiosTable;


    // This is the data that is stored and restored.

    // coalCounts[b][j] = number of coalescences in branch b in gtree j = k_jb in paper
    // coalIntensities[b][j] = gamma_bj = (1/p_j) sum_{i=0}^{k_jb} c_jbi (n_jb-i choose 2)
    // sum_j gamma_bj = gamma_b in paper.
    private int [][] coalCounts;
    private double [][] coalIntensities;

    // The per-gene info in bindings is only fully updated when all gene trees are compatible,
    // and the compatibility is not recalculated beyond finding one incompatible gene tree.
    // The two arrays of dirty flags are used to deal with this.
    // gTreeFitIsDirty[] is per-gene. Says whether we know if a gtree fits stree.
    // gTreeFits[] is per-gene. If (gTreeFitIsDirty[j] == false), this says whether a gtree fits stree.
    // Furthermore, if (gTreeFitIsDirty[j] == false & gTreeFits[j] == true), it means the heights
    // of coalescences in each stree branch have been found.
    private boolean [] gTreeFitIsDirty;
    private boolean [] gTreeFits;

    // gTreeCountIntensityIsDirty[j] says whether coalCounts[][j] and
    // coalIntensities[][j] are up to date for gtree j (for all stree nodes).
    private boolean [] gTreeCountIntensityIsDirty;


    // These are the stored version of the above.
    private int [][] storedCoalCounts;
    private double [][] storedCoalIntensities;

    private boolean [] storedGTreeFits;
    private boolean [] storedGTreeFitIsDirty;

    private boolean [] storedGTreeCountIntensityIsDirty;


    private boolean debugFlag = Boolean.valueOf(System.getProperty("stacey.debug"));
    private int numberofdebugchecks = 0;
    private final static int maxnumberofdebugchecks = 1000000;



    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate(); // ??

        gTreeCFs = geneTreesInput.get();
        nGTrees = gTreeCFs.size();

        gTrees = new ArrayList<>(nGTrees);
        for (int j = 0; j < nGTrees; j++) {
            gTrees.add(gTreeCFs.get(j).getTree());
        }

        sTree = treeInput.get();
        bindings = Bindings.initialise(sTree, gTrees);
        fitsHeights = new FitsHeights(sTree, gTrees, bindings);
        nSMCTreeNodes = sTree.getNodeCount();
        nSMCTreeTips = sTree.getLeafNodeCount();

        int totalGTipCount = 0;
        for (Tree gtree : gTrees) {
            totalGTipCount += gtree.getLeafNodeCount();
        }

        lnGammaRatiosTable = new double[priorComponentsInput.get().size()][totalGTipCount];
        for (int c = 0; c < priorComponentsInput.get().size(); c++) {
            InverseGammaComponent igc = priorComponentsInput.get().get(c);
            for (int q = 0; q < totalGTipCount; q++) {
                lnGammaRatiosTable[c][q] = lnRatioGammas(igc.getAlpha(), q);
            }
        }

        coalCounts = new int[nSMCTreeNodes][nGTrees];
        coalIntensities = new double[nSMCTreeNodes][nGTrees];
        storedCoalCounts = new int[nSMCTreeNodes][nGTrees];
        storedCoalIntensities = new double[nSMCTreeNodes][nGTrees];

        // normalize the inverse gamma component weights to sum to 1
        double totalweight = 0.0;
        for (int i = 0; i < priorComponentsInput.get().size(); i++) {
            totalweight += priorComponentsInput.get().get(i).getWeight();
        }
        for (int i = 0; i < priorComponentsInput.get().size(); i++) {
            priorComponentsInput.get().get(i).normalizeWeight(totalweight);
        }
        // set up all the dirtiness flags
        gTreeFits =    new boolean[nGTrees];
        storedGTreeFits = new boolean[nGTrees];

        gTreeFitIsDirty =    new boolean[nGTrees];
        storedGTreeFitIsDirty = new boolean[nGTrees];

        for (int j = 0; j < nGTrees; j++) {
            gTreeFitIsDirty[j] =    true;
            storedGTreeFitIsDirty[j] = true;
        }

        gTreeCountIntensityIsDirty = new boolean [nGTrees];
        storedGTreeCountIntensityIsDirty = new boolean [nGTrees];
        for (int j = 0; j < nGTrees; j++) {
            gTreeCountIntensityIsDirty[j] =    true;
            storedGTreeCountIntensityIsDirty[j] = true;
        }
    }



    @Override
    public double calculateLogP() {
        // Basically, this calls logLhoodAllGeneTreesInSMCTree(); the rest is debug checks

        double fastLogP = logLhoodAllGeneTreesInSMCTree(getInverseGammaMixture(),
                popPriorScale.get().getValue(), false);

        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            double robustLogP = logLhoodAllGeneTreesInSMCTree(getInverseGammaMixture(),
                    popPriorScale.get().getValue(), true);
            if (Math.abs(fastLogP - robustLogP) > 1e-12) {
                System.out.println("BUG in calculateLogP() in PIOMSCoalescentDistribution");
                assert false;
            }
            numberofdebugchecks++;
        }

        if (debugFlag && numberofdebugchecks < 3) {
            System.out.println("debugging: coal logP " + fastLogP);
            System.out.println(Misc.allTreesAsText(sTree, gTrees.subList(0,2)));
        }
        logP = fastLogP;
        return logP;
    }



    @Override
    protected boolean requiresRecalculation() {
        if (sTree.somethingIsDirty()) {
            for (int j = 0; j < nGTrees; j++) {
                gTreeFitIsDirty[j] = gTreeCountIntensityIsDirty[j] = true;
            }
        }
        for (int j = 0; j < nGTrees; j++) {
            boolean gTreeDirty = gTrees.get(j).somethingIsDirty();
            // I don't update everything in case of incompatibility, so need to OR
            gTreeFitIsDirty[j] |= gTreeDirty;
            gTreeCountIntensityIsDirty[j] |=  gTreeDirty;
        }
        return true;
    }


    @Override
    public List<String> getArguments() {
        return null;
    }


    @Override
    public List<String> getConditions() {
        return null;
    }


    @Override
    public void sample(final State state, final Random random) {
    }


    @Override
    public void store() {
        System.arraycopy(gTreeFits, 0, storedGTreeFits, 0, storedGTreeFits.length);
        System.arraycopy(gTreeFitIsDirty, 0, storedGTreeFitIsDirty, 0, storedGTreeFitIsDirty.length);
        for (int n = 0; n < nSMCTreeNodes; n++) {
            System.arraycopy(coalCounts[n], 0, storedCoalCounts[n], 0, storedCoalCounts[n].length);
            System.arraycopy(coalIntensities[n], 0, storedCoalIntensities[n], 0, storedCoalIntensities[n].length);
        }
        System.arraycopy(gTreeCountIntensityIsDirty, 0, storedGTreeCountIntensityIsDirty, 0, storedGTreeCountIntensityIsDirty.length);
        super.store();
    }



    @Override
    public void restore() {
        System.arraycopy(storedGTreeFits, 0, gTreeFits, 0, gTreeFits.length);
        int [][] tmpCCs = coalCounts;
        coalCounts = storedCoalCounts;
        storedCoalCounts = tmpCCs;
        double [][] tmpCIs = coalIntensities;
        coalIntensities = storedCoalIntensities;
        storedCoalIntensities = tmpCIs;
        System.arraycopy(storedGTreeCountIntensityIsDirty, 0, gTreeCountIntensityIsDirty, 0, gTreeCountIntensityIsDirty.length);
        super.restore();
    }


    /***************************************************************************************************/
    /***************************************************************************************************/
    /***************************************************************************************************/


    // used here and by PopSampleStatistic
    InverseGammaMixture getInverseGammaMixture() {
        int n = priorComponentsInput.get().size();
        InverseGammaMixture igm = new InverseGammaMixture(n);
        for (int c = 0; c < n; c++) {
            igm.setWeight(c, priorComponentsInput.get().get(c).getWeight());
            igm.setAlpha(c, priorComponentsInput.get().get(c).getAlpha());
            igm.setBeta(c, priorComponentsInput.get().get(c).getBeta());
        }
        return igm;
    }


    /***************************************************************************************************/
    /*******************************           private              ************************************/
    /***************************************************************************************************/

    private double logLhoodAllGeneTreesInSMCTree(InverseGammaMixture igm,
                                                 double popPriorScale, boolean robust) {

        if (robust) {
            for (int j = 0; j < nGTrees; j++) {
                gTreeFitIsDirty[j] = gTreeCountIntensityIsDirty[j] = true;
            }
        }

        /*
        if (debugFlag   &&  numberofdebugchecks < maxnumberofdebugchecks &&  !robust) {
            if (numberofdebugchecks < 10  ||  numberofdebugchecks == 52 ||
                    (numberofdebugchecks < 100  &&  numberofdebugchecks % 10 == 0)  ||
                    (numberofdebugchecks < 1000  &&  numberofdebugchecks % 100 == 0) ||
                    (numberofdebugchecks < 10000  &&  numberofdebugchecks % 1000 == 0) ||
                    (numberofdebugchecks % 100000 == 0)) {
                System.out.println("numberofdebugchecks = " + numberofdebugchecks);
                System.out.println(Misc.beastSubtreeAsTextWithHeader(sTree.getRoot(), null));
                for (int j = 0; j < gTrees.size(); j++) {
                    System.out.println(Misc.beastSubtreeAsTextWithHeader(gTrees.get(j).getRoot(), null));
                }
                System.out.println("");
                System.out.flush();
            }
        }*/
        
        // unless gtree fit is clean and the gtree fits, CountIntensity dirtyness may be wrong
        for (int j = 0; j < nGTrees; j++) {
            if (gTreeFitIsDirty[j]  ||  !gTreeFits[j]) {
                gTreeCountIntensityIsDirty[j] = true;
            }
        }
        // make gTreeFits up to date
        boolean allGtreesFit = true;
        for (int j = 0; j < nGTrees  &&  allGtreesFit; j++) {
            if (gTreeFitIsDirty[j]) {
                gTreeFits[j] = fitsHeights.updateFitHeightsForOneGTree(j);
                allGtreesFit = allGtreesFit && gTreeFits[j];
                gTreeFitIsDirty[j] = false;
            } else {
                if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
                    // consistency check
                    if (gTreeFits[j] != fitsHeights.updateFitHeightsForOneGTree(j)) {
                        System.out.println("BUG in logLhoodAllGeneTreesInSMCTree() gTreeFits[j] wrong");
                        System.exit(1);
                        assert false;
                    }
                }
            }
        }
        // if incompatibility, return -oo
        if (!allGtreesFit) {
            return Double.NEGATIVE_INFINITY;
        }
        // find max height of any gene tree for intensity at stree root
        double maxgtreehgt = 0.0;
        for (Tree gTree : gTrees) {
            maxgtreehgt = Math.max(maxgtreehgt, gTree.getRoot().getHeight());
        }
        // make Nlineages and coal counts and intensities up to date
        for (int j = 0; j < nGTrees; j++) {
            if (gTreeCountIntensityIsDirty[j]) {
                // counts
                double coalFactor = gTreeCFs.get(j).getCoalFactor();
                for (int n = 0; n < nSMCTreeNodes; n++) {
                    coalCounts[n][j] = fitsHeights.getHeightsFromSNodeNrGTree(n, j).size();
                }
                // nlineages
                int [] nLins = new int[nSMCTreeNodes];
                for (int n = 0; n < nSMCTreeTips; n++) {
                    nLins[n] = bindings.nLineagesForBeastTipNrAndGtree(n, j);
                }
                fillinSubtreeNLineages(nLins, sTree.getRoot(), j);
                // intensities
                for (int n = 0; n < nSMCTreeNodes; n++) {
                    ArrayList<Double> njHeights = fitsHeights.getHeightsFromSNodeNrGTree(n, j);

                    double nodeHeight = sTree.getNode(n).getHeight();
                    if (debugFlag) {
                        for (int i = 0; i < njHeights.size(); i++) {
                            if (njHeights.get(i) < nodeHeight) {
                                System.out.println("BUG in logLhoodAllGeneTreesInSMCTree() bad height (1).");
                            }
                        }
                    }

                    Node anc = sTree.getNode(n).getParent();
                    double ancHeight = (anc == null) ? maxgtreehgt : anc.getHeight();
                    int nCoals = njHeights.size();

                    assert nCoals == coalCounts[n][j];
                    assert nCoals < nLins[n];
                    if (debugFlag) {
                        for (int i = 0; i < njHeights.size(); i++) {
                            if (njHeights.get(i) > ancHeight) {
                                System.out.println("BUG in logLhoodAllGeneTreesInSMCTree() bad height (2).");
                            }
                        }
                    }

                    double [] heights = new double[nCoals+2];
                    heights[0] = nodeHeight;
                    heights[nCoals+1] = ancHeight;
                    for (int i = 0; i < nCoals; i++) {
                        heights[i+1] = njHeights.get(i);
                    }
                    if (nCoals >= 2) {
                        Arrays.sort(heights);
                    }
                    double coalIntensity = 0.0; // gamma_nj
                    for (int i = 0; i < nCoals+1; i++) {
                        coalIntensity += (heights[i+1] - heights[i]) * (nLins[n] - i) * (nLins[n] - i - 1) * 0.5;
                    }
                    coalIntensity /= coalFactor;

                    // old calc
/*                    if (nCoals == 0) {
                        coalIntensities[n][j] = (ancHeight - nodeHeight) * nLins[n] * (nLins[n] - 1) * 0.5 / coalFactor;
                    } else if (nCoals == 1) {
                        double gammanj = (njHeights.get(0) - nodeHeight) * nLins[n] * (nLins[n] - 1) * 0.5;
                        gammanj += (ancHeight - njHeights.get(nCoals-1)) * (nLins[n] - nCoals) * (nLins[n] - nCoals - 1) * 0.5;
                        coalIntensities[n][j] = gammanj / coalFactor;
                    } else {
                        njHeights.sort(null);
                        double gammanj = (njHeights.get(0) - nodeHeight) * nLins[n] * (nLins[n] - 1) * 0.5;
                        for (int h = 1; h < nCoals; h++) {
                            gammanj += (njHeights.get(h) - njHeights.get(h-1)) * (nLins[n] - h) * (nLins[n] - h - 1) * 0.5;
                        }
                        gammanj += (ancHeight - njHeights.get(nCoals-1)) * (nLins[n] - nCoals) * (nLins[n] - nCoals - 1) * 0.5;
                        coalIntensities[n][j] = gammanj / coalFactor;
                    }
                    // check
                    if (coalIntensities[n][j] != coalIntensity) {
                        System.out.println("BUG in logLhoodAllGeneTreesInSMCTree() check coalIntensity failed");
                    }
                    assert coalIntensities[n][j] == coalIntensity;*/
                    coalIntensities[n][j] = coalIntensity;
                }
                gTreeCountIntensityIsDirty[j] = false;
            }
        }
        // calculate and sanity-check the result
        double logPGS = logProbAllGTreesInSMCTree(igm, popPriorScale);
        assert !Double.isNaN(logPGS);
        assert !Double.isInfinite(logPGS);
        return logPGS;
    }


    private void fillinSubtreeNLineages(int [] nlineages, Node node, int j) {
        int lftNLin;
        Node lftNode = node.getChild(0);
        if (lftNode.isLeaf()) {
            lftNLin = nlineages[lftNode.getNr()];
        } else {
            fillinSubtreeNLineages(nlineages, lftNode, j);
            lftNLin = nlineages[lftNode.getNr()];
        }
        lftNLin -= coalCounts[lftNode.getNr()][j];
        int rgtNLin;
        Node rgtNode = node.getChild(1);
        if (rgtNode.isLeaf()) {
            rgtNLin = nlineages[rgtNode.getNr()];
        } else {
            fillinSubtreeNLineages(nlineages, rgtNode, j);
            rgtNLin = nlineages[rgtNode.getNr()];
        }
        rgtNLin -= coalCounts[rgtNode.getNr()][j];
        nlineages[node.getNr()] = lftNLin + rgtNLin;
    }



    private double logProbAllGTreesInSMCTree(InverseGammaMixture igm, double popPriorScale) {
        double logP = 0.0;
        
        double [] logCFs = new double[nGTrees];
        for (int j =0; j < nGTrees; j++) {
            logCFs[j] = Math.log(gTreeCFs.get(j).getCoalFactor());
        }
        
        for (int n = 0; n < nSMCTreeNodes; n++) {
            int q_b = 0;
            double gamma_b = 0.0;
            double minusLog_r_b = 0.0;
            for (int j= 0; j < nGTrees; j++) {
                q_b += coalCounts[n][j];
                minusLog_r_b += logCFs[j] * coalCounts[n][j];
                gamma_b += coalIntensities[n][j];
            }
            double [] lambdas = igm.getWeights();
            double [] alphas = igm.getAlphas();
            double [] betas = igm.getBetas();
            double [] logProbCpts = new double[lambdas.length];

            for (int c = 0; c < lambdas.length; c++) {
                logProbCpts[c] = 0.0;
                double sigmabeta_c = popPriorScale * betas[c];
                logProbCpts[c] += Math.log(lambdas[c]);
                logProbCpts[c] += alphas[c] * Math.log(sigmabeta_c);
                logProbCpts[c] -= (alphas[c] + q_b) * Math.log(sigmabeta_c + gamma_b);
                logProbCpts[c] += lnGammaRatiosTable[c][q_b];
            }
            double logP_b = logSumExp(logProbCpts);
            logP_b -= minusLog_r_b;
            logP += logP_b;
        }
        return logP;
    }



    // returns log(Gamma(a+n)/Gamma(a))
    private double lnRatioGammas(double a, int n) {
        double x = 0.0;
        for (int i = 0; i < n; i++) {
            x += Math.log(a + i);
        }
        return x;
    }


    // Returns log(sum_i e^(x_i)) without over/under flow
    private double logSumExp(double x[]) {
        double maxx = Double.NEGATIVE_INFINITY;
        for (double d : x) {
            if (d > maxx) { maxx = d; }
        }
        double sum = 0.0;
        for (double d : x) {
            sum += Math.exp(d-maxx);
        }
        return maxx + Math.log(sum);
    }

}

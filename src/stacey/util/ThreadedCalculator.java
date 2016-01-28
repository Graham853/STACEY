package stacey.util;

        import beast.evolution.tree.Node;
        import beast.evolution.tree.Tree;
        import beast.evolution.tree.TreeInterface;
        import stacey.GtreeAndCoalFactor;

        import java.util.ArrayList;
        import java.util.Arrays;
        import java.util.List;
        import java.util.concurrent.CountDownLatch;
        import java.util.concurrent.ExecutorService;
        import java.util.concurrent.Executors;
        import java.util.concurrent.atomic.AtomicBoolean;

/**
 * Created by Work on 21/12/2015.
 */
public class ThreadedCalculator {
    private ExecutorService exec = null;
    private CountDownLatch countDown;
    private static AtomicBoolean allGtreesFit = new AtomicBoolean(true);

    public ThreadedCalculator(int nThreads) {
        exec = Executors.newFixedThreadPool(nThreads);
    }


    // The runnable class.
    class GTreeLhoodInfoIntoSTreeRunner implements Runnable {
        PIOMSCInfoForLogPCalc info;
        int j;

        GTreeLhoodInfoIntoSTreeRunner(PIOMSCInfoForLogPCalc info, int j) {
            this.info = info;
            this.j = j;
        }

        public void run() {
            try {
                putOneGTreeLhoodInfoIntoSTree(info, j);
            } catch (Exception e) {
                System.err.println("Error in GTreeLhoodInfoIntoSTreeRunner().");
                throw new RuntimeException(e.getMessage());
            }
            countDown.countDown();
        }
    }



    public void putAllGTreeLhoodInfosIntoSTree(PIOMSCInfoForLogPCalc info) {
        int nGTrees = info.getGTreeFits().length;
        try {
            countDown = new CountDownLatch(nGTrees);
            for (int j = 0; j < nGTrees; j++) {
                GTreeLhoodInfoIntoSTreeRunner runner = new GTreeLhoodInfoIntoSTreeRunner(info, j);
                exec.execute(runner);
            }
            countDown.await();
        } catch (Exception e) {
            System.err.println("Error in putAllGTreeLhoodInfosIntoSTree().");
            throw new RuntimeException(e.getMessage());
        }
    }


    public boolean getAllGtreesFit() {
        return allGtreesFit.get();
    }




    private static void putOneGTreeLhoodInfoIntoSTree(PIOMSCInfoForLogPCalc info, int j) {
        FitsHeights fitsHeights = info.getFitsHeights(); // update and use the j'th part
        Bindings bindings = info.getBindings(); // used for nLineages at tips
        TreeInterface sTree = info.getSTree(); // used for nLineages, node heights
        List<GtreeAndCoalFactor> gTreeCFs = info.getGTreeCFs(); // used
        double maxgtreehgt = info.getMaxgtreehgt(); // used
        boolean [] gTreeFitIsDirty = info.getGTreeFitIsDirty(); // [j]. update and use the j'th part
        boolean [] gTreeFits = info.getGTreeFits(); // [j]. update and use the j'th part
        boolean [] gTreeCountIntensityIsDirty = info.getGTreeCountIntensityIsDirty(); // [j]. update and use the j'th part
        int [][] coalCounts = info.getCoalCounts(); // [n][j]. update and use the j'th part for all n
        double [][] coalIntensities = info.getCoalIntensities(); // [n][j]. update and use the j'th part for all n
        int [] q_j = info.getQ_j(); // [n]. update
        double [] gamma_j = info.getGamma_j(); // [n]. update
        double [] minusLog_r_j = info.getMinusLog_r_j(); // [n]. update

        // check that no other thread has found incompatibility: if it has, return
        if (!allGtreesFit.get()) {
            return;
        }
        // update gTreeFits[j]. If it's false, set allGtreesFit and return.
        if (gTreeFitIsDirty[j]) {
            gTreeFits[j] = fitsHeights.updateFitHeightsForOneGTree(j);
            if (!gTreeFits[j]) {
                allGtreesFit.compareAndSet(true, false);
            }
            gTreeFitIsDirty[j] = false;
            if (!gTreeFits[j]) {
                return;
            }
        }
        // gTree j fits, so continue
        int nSMCTreeNodes = sTree.getNodeCount();
        int nSMCTreeTips = sTree.getLeafNodeCount();
        // Make coalCounts, coalIntensities up to date
        if (gTreeCountIntensityIsDirty[j]) {
            // counts
            for (int b = 0; b < nSMCTreeNodes; b++) {
                coalCounts[b][j] = fitsHeights.getHeightsFromSNodeNrGTree(b, j).size();
            }
            // nlineages
            int [] nLins = new int[nSMCTreeNodes];
            for (int n = 0; n < nSMCTreeTips; n++) {
                nLins[n] = bindings.nLineagesForBeastTipNrAndGtree(n, j);
            }
            fillinSubtreeNLineages(nLins, coalCounts, sTree.getRoot(), j);
            // intensities
            for (int b = 0; b < nSMCTreeNodes; b++) {
                ArrayList<Double> njHeights = fitsHeights.getHeightsFromSNodeNrGTree(b, j);
                double nodeHeight = sTree.getNode(b).getHeight();
                Node anc = sTree.getNode(b).getParent();
                double ancHeight = (anc == null) ? maxgtreehgt : anc.getHeight();
                int nCoals = njHeights.size();
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
                    coalIntensity += (heights[i+1] - heights[i]) * (nLins[b] - i) * (nLins[b] - i - 1) * 0.5;
                }
                coalIntensity /= gTreeCFs.get(j).getCoalFactor();
                coalIntensities[b][j] = coalIntensity;
            }
            gTreeCountIntensityIsDirty[j] = false;
        }
        // Fill in q_j, gamma_j, minusLog_r_j.
        for (int b = 0; b < nSMCTreeNodes; b++) {
            q_j[b] += coalCounts[b][j];
            minusLog_r_j[b] += Math.log(gTreeCFs.get(j).getCoalFactor()) * coalCounts[b][j];
            gamma_j[b] += coalIntensities[b][j];
        }
    }



    // called from putOneGTreeLhoodInfoIntoSTree()
    private static void fillinSubtreeNLineages(int [] nlineages, int [][] coalCounts, Node node, int j) {
        int lftNLin;
        Node lftNode = node.getChild(0);
        if (lftNode.isLeaf()) {
            lftNLin = nlineages[lftNode.getNr()];
        } else {
            fillinSubtreeNLineages(nlineages, coalCounts, lftNode, j);
            lftNLin = nlineages[lftNode.getNr()];
        }
        lftNLin -= coalCounts[lftNode.getNr()][j];
        int rgtNLin;
        Node rgtNode = node.getChild(1);
        if (rgtNode.isLeaf()) {
            rgtNLin = nlineages[rgtNode.getNr()];
        } else {
            fillinSubtreeNLineages(nlineages, coalCounts, rgtNode, j);
            rgtNLin = nlineages[rgtNode.getNr()];
        }
        rgtNLin -= coalCounts[rgtNode.getNr()][j];
        nlineages[node.getNr()] = lftNLin + rgtNLin;
    }





    // TODO this to go in PIOMSCoalescentDistribution. It is entirely untested
    /*
    private double ThreadedLogP(InverseGammaMixture igm,
                                double popPriorScale) {
        ThreadedCalculator thcc = new ThreadedCalculator(4); // TODO call in initAndValidate() ?

        // find max height of any gene tree for intensity at stree root
        double maxgtreehgt = 0.0;
        for (Tree gTree : gTrees) {
            maxgtreehgt = Math.max(maxgtreehgt, gTree.getRoot().getHeight());
        }

        int [] q_j = new int[sTree.getNodeCount()];
        double [] gamma_j = new double[sTree.getNodeCount()];
        double [] minusLog_r_j = new double[sTree.getNodeCount()];
        PIOMSCInfoForLogPCalc info =
                new PIOMSCInfoForLogPCalc(fitsHeights, bindings, sTree, gTreeCFs, maxgtreehgt,
                        gTreeFitIsDirty, gTreeFits, gTreeCountIntensityIsDirty,
                        coalCounts, coalIntensities,
                        q_j,  gamma_j, minusLog_r_j);
        thcc.putAllGTreeLhoodInfosIntoSTree(info);
        // threads done
        // if incompatibility, return -oo
        if (!thcc.getAllGtreesFit()) {
            return Double.NEGATIVE_INFINITY;
        }
        double [] lambdas = igm.getWeights();
        double [] alphas = igm.getAlphas();
        double [] betas = igm.getBetas();
        double [] logProbCpts = new double[lambdas.length];
        double logPGS = 0.0;
        for (int b = 0; b < sTree.getNodeCount(); b++) {
            for (int c = 0; c < lambdas.length; c++) {
                logProbCpts[c] = 0.0;
                double sigmabeta_c = popPriorScale * betas[c];
                logProbCpts[c] += Math.log(lambdas[c]);
                logProbCpts[c] += alphas[c] * Math.log(sigmabeta_c);
                logProbCpts[c] -= (alphas[c] + q_j[b]) * Math.log(sigmabeta_c + gamma_j[b]);
                logProbCpts[c] += lnGammaRatiosTable[c][q_j[b]];
            }
            double logP_b = logSumExp(logProbCpts);
            logP_b -= minusLog_r_j[b];
            logPGS += logP_b;
        }

        assert !Double.isNaN(logPGS);
        assert !Double.isInfinite(logPGS);
        return logPGS;
    }*/





}

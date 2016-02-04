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

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;
import stacey.debugtune.Checks;
import stacey.util.*;

import java.util.ArrayList;
import java.util.List;


@Description("A move which scales some node heights in the SMC-tree and gene trees. " +
        "The further away a node is from a 'focal' SMC-tree node, the less it is affected.")
public class FocusedNodeHeightScaler extends Operator {

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<Tree> smcTreeInput =
            new Input<Tree>("smcTree",
                    "The species tree or minimal clusters tree", Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<List<Tree>> geneTreesInput =
            new Input<List<Tree>>("geneTree",
                    "All gene trees",
                    new ArrayList<Tree>());

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<Long> delayInput =
            new Input<>("delay",
                    "Number of times the operator is disabled.");

    private Tree sTree;
    private List<Tree> gTrees;
    private long delay = 0;
    private int callCount = 0;
    private boolean sTreeTooSmall;

    private UnionArrays unionArrays;

    private double [] sTreeLogHeights;
    private int [] sTreeDistances;
    private double [] sTreeWeights;

    private double [][] gTreeLogHeights;
    private int [][] gTreeDistances;
    private double [][] gTreeWeights;


    private final boolean debugFlag = Boolean.valueOf(System.getProperty("stacey.debug"));
    private int numberofdebugchecks = 0;
    private final static int maxnumberofdebugchecks = 100000;




    /* This is used to pass information from the smcTree node chosen as
    focus to the gene trees so they can calculate distances from the focus.
    This is used first, FSMoveSMCNodeInfo is used later in the calculation
    for the move. */
    private static class OpFSinfoSMCNodeUnions {
        private final BitUnion nodeUnion;
        private final BitUnion lftUnion;
        private final BitUnion rgtUnion;

        public OpFSinfoSMCNodeUnions(BitUnion nodeUnion, BitUnion lftUnion, BitUnion rgtUnion) {
            this.nodeUnion = nodeUnion;
            this.lftUnion = lftUnion;
            this.rgtUnion = rgtUnion;

        }

        public BitUnion getNodeUnion() { return nodeUnion; }
        public BitUnion getLftUnion() { return lftUnion; }
        public BitUnion getRgtUnion() { return rgtUnion; }

    }



    /* This is used to pass information from all smcTree nodes to the gene trees
    so they can calculate bounds for the scaling. */
    private static class OpFSinfoSMCNode {
        private final BitUnion lftUnion;
        private final BitUnion rgtUnion;
        private final double logHeight;
        private final double weight;

        public OpFSinfoSMCNode(BitUnion nodeUnion, BitUnion lftUnion, BitUnion rgtUnion,
                               double logHeight, double weight) {
            this.lftUnion = lftUnion;
            this.rgtUnion = rgtUnion;
            this.logHeight = logHeight;
            this.weight = weight;
        }

        public BitUnion getLftUnion() { return lftUnion; }
        public BitUnion getRgtUnion() { return rgtUnion; }
        public double getLogHeight() { return logHeight; }
        public double getWeight() { return weight; }
    }


    /*******************************************************************************/


    @Override
    public void initAndValidate() {
        /*if (smcTreeInput.get().getLeafNodeCount() < 5) {
            throw new Exception("FocusedNodeHeightScaler cannot be used if there are less than 4 minimal clusters.");
        } TODO this upsets Beauti */

        sTree = smcTreeInput.get();
        gTrees = geneTreesInput.get();
        sTreeTooSmall = (sTree.getLeafNodeCount() < 5);
        if (delayInput != null  &&  delayInput.get() != null) {
            delay = delayInput.get().longValue();
        }
        Bindings bindings = Bindings.initialise(sTree, gTrees);

        unionArrays = UnionArrays.initialise(sTree, gTrees, bindings);

        sTreeLogHeights = new double[smcTreeInput.get().getNodeCount()];
        sTreeDistances = new int[smcTreeInput.get().getNodeCount()];
        sTreeWeights = new double[smcTreeInput.get().getNodeCount()];

        gTreeLogHeights = new double[gTrees.size()][];
        gTreeDistances = new int[gTrees.size()][];
        gTreeWeights = new double[gTrees.size()][];
        for (int j = 0; j < gTrees.size(); j++) {
            TreeInterface gTree = gTrees.get(j);
            gTreeLogHeights[j] = new double[gTree.getNodeCount()];
            gTreeDistances[j] = new int[gTree.getNodeCount()];
            gTreeWeights[j] = new double[gTree.getNodeCount()];
        }
    }



    @Override
    public double proposal() {
        // TODO Input says use get(this). See TODO in CoordinatedPruneRegraft
        if (sTreeTooSmall) {
            return Double.NEGATIVE_INFINITY;
        }
        callCount++;
        if (callCount < delay) {
            return Double.NEGATIVE_INFINITY;
        }


        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            Checks.allTreesAndCompatibility(sTree, gTrees, "FocusedNodeHeightScaler", "before move");
            numberofdebugchecks++;
        }

        unionArrays.update();
        double logHR = doFocusedScalerMove();   // The business
        unionArrays.reset();

        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            Checks.allTreesAndCompatibility(sTree, gTrees, "FocusedNodeHeightScaler", "after move");
            numberofdebugchecks++;
        }
        return logHR;
    }

    /****************************************************************************************/


    private double doFocusedScalerMove() {

        sTree.startEditing(this);
        // TODO This is a fix for beast.core.State$Trie memory
        for (int j = 0; j < gTrees.size(); j++) {
            gTrees.get(j).startEditing(this);
        }

        // Stage 1 - choose focal node and fill in info.
        // Fill in log(heights) and distances for smc tree nodes.
        // Choose focal node, get unions of itself and children.
        // Convert distances to weights for the smc tree.
        // Fill in log(heights) and distances for each gene tree.
        // Convert distances to weights for each gene tree's nodes
        OpFSinfoSMCNodeUnions fsNUnions = setUpForFSMove(sTreeLogHeights, sTreeDistances);
        int rootDistance = sTreeDistances[sTreeDistances.length-1];
        for (int n = 0; n < sTreeDistances.length; n++) {
            sTreeWeights[n] = distanceToWeight(sTreeDistances[n], rootDistance);
        }
        // TODO-threaded
        for (int j = 0; j < gTrees.size(); j++) {
            int [] jDistances =  gTreeDistances[j];
            double [] jLogHeights = gTreeLogHeights[j];
            fillInGTreeDistancesLogHeights(j, fsNUnions, jDistances, jLogHeights);
            int jRootDistance = jDistances[jDistances.length-1];
            double [] jWeights = gTreeWeights[j];
            for (int n = 0; n < jDistances.length; n++) {
                jWeights[n] = distanceToWeight(jDistances[n], jRootDistance);
            }
        }

        // Stage 2 - find bounds.
        // Bounds from branches within smc tree.
        // Bounds from branches within gene trees and between smc tree and gene trees.
        double [] interval = new double[2];
        interval[0] = Double.NEGATIVE_INFINITY;
        interval[1] = Double.POSITIVE_INFINITY;
        getSTreeFSBounds(interval, sTreeLogHeights, sTreeWeights);
        OpFSinfoSMCNode[] fsNodeInfos = getAllInternalFSMoveSNInfos(sTreeLogHeights, sTreeWeights);
        // TODO-threaded
        for (int j = 0; j < gTrees.size(); j++) {
            double [] jLogHeights = gTreeLogHeights[j];
            double [] jWeights = gTreeWeights[j];
            getGTreeFSBounds(interval, j, fsNodeInfos, jLogHeights, jWeights);
        }
        assert !Double.isInfinite(interval[0]);
        assert !Double.isInfinite(interval[1]);
        assert interval[0] < 0.0;
        assert interval[1] > 0.0;

        // Stage 3 - carry out the move
        double logHR = 0.0;
        double nearlyOne = 0.99999999;  // to avoid numerical errors which might make a branch slightly negative
        double logSF = Randomizer.uniform(nearlyOne*interval[0], nearlyOne*interval[1]);
        logHR += fsMoveDoScale(sTree, logSF, sTreeWeights);
        // TODO-threaded
        for (int j = 0; j < gTrees.size(); j++) {
            double [] jWeights = gTreeWeights[j];
            logHR += fsMoveDoScale(gTrees.get(j), logSF, jWeights);
        }
        return logHR;
    }



    /*************************************************************************************/
    /*************************   Dealing with SMC tree ***********************************/
    /*************************************************************************************/


    private OpFSinfoSMCNodeUnions setUpForFSMove(double [] logHeights, int [] distances) {
        fillInAnyTreeLogHeights(sTree, logHeights);
        int s;
        do {
            s = Randomizer.nextInt(sTree.getNodeCount());
        } while (!nodeCanBeFocus(s));
        Node sN = sTree.getNode(s);
        BitUnion sunion = unionArrays.sNodeUnion(s);
        BitUnion lftunion = unionArrays.sNodeUnion(sN.getChild(0).getNr());
        BitUnion rgtunion = unionArrays.sNodeUnion(sN.getChild(1).getNr());
        fillInSTreeDistances(sN, distances);
        return new OpFSinfoSMCNodeUnions(sunion, lftunion, rgtunion);
    }


    private boolean nodeCanBeFocus(int s) {
        // not a tip, not the root, and at least one child not a tip
        if (sTree.getNode(s).isLeaf()) {
            return false;
        }if (sTree.getNode(s).isRoot()) {
            return false;
        }
        Node lftN = sTree.getNode(s).getChild(0);
        Node rgtN = sTree.getNode(s).getChild(1);
        return (!lftN.isLeaf() ||  !rgtN.isLeaf());
    }



    private void getSTreeFSBounds(double[] interval, double[] logHeights, double[] weights) {
        getAnyTreeInternalFSBounds(interval, sTree, logHeights, weights);
    }


    private void fillInSTreeDistances(Node node, int[] distances) {
        assert sTree.getNodeCount() == distances.length;
        for (int i = 0; i < distances.length; i++) {
            distances[i] = Integer.MAX_VALUE;
        }
        distances[node.getNr()] = 0;
        while (!node.isRoot()) {
            Node ancN = node.getParent();
            int anc = ancN.getNr();
            distances[anc] = Math.min(distances[anc], distances[node.getNr()] + 1);
            node = ancN;
        }
        fillInAnySubtreeDistances(sTree.getRoot(), distances);
    }


    private OpFSinfoSMCNode[] getAllInternalFSMoveSNInfos(double[] logHeights, double[] weights) {
        OpFSinfoSMCNode[] fsNInfos = new OpFSinfoSMCNode[sTree.getInternalNodeCount()];
        int i = 0;
        for (int s = 0; s < sTree.getNodeCount(); s++) {
            Node node = sTree.getNode(s);
            if (!node.isLeaf()) {
                assert i < fsNInfos.length;
                Node lftN = node.getChild(0);
                Node rgtN = node.getChild(1);
                BitUnion sunion = unionArrays.sNodeUnion(s);
                BitUnion lftunion = unionArrays.sNodeUnion(lftN.getNr());
                BitUnion rgtunion = unionArrays.sNodeUnion(rgtN.getNr());
                double logHeight = logHeights[s];
                double weight = weights[s];
                fsNInfos[i++] = new OpFSinfoSMCNode(sunion, lftunion, rgtunion, logHeight, weight);
            }
        }
        assert i == fsNInfos.length;
        return fsNInfos;
    }




    /*************************************************************************************/
    /*************************   Dealing with gene trees *********************************/
    /*************************************************************************************/



    private void fillInGTreeDistancesLogHeights(int j, OpFSinfoSMCNodeUnions fsNUnions, int[] distances, double[] logHeights) {
        fillInGTreeDistances(j, fsNUnions, distances);
        fillInGTreeLogHeights(j, logHeights);
    }



    private void getGTreeFSBounds(double[] interval, int j, OpFSinfoSMCNode[] fsNodeInfos,
                                 double[] logHeights, double[] weights) {
        getGTreeFocusedScalerInternalBounds(interval, j, logHeights, weights);
        for (int i = 0; i < fsNodeInfos.length; i++) {
            getGTreeFocusedScalerCompatibilityBounds(interval, j, fsNodeInfos[i], logHeights, weights);
        }
    }




    private void fillInGTreeDistances(int j, OpFSinfoSMCNodeUnions fsNUnions, int [] distances) {
        TreeInterface gTree = gTrees.get(j);
        assert gTree.getNodeCount() == distances.length;
        BitUnion spp = fsNUnions.getNodeUnion();
        BitUnion sppl = fsNUnions.getLftUnion();
        BitUnion sppr = fsNUnions.getRgtUnion();
        for (int i = 0; i < gTree.getNodeCount(); i++) {
            if (gNodeSNodeCloselyLinked(j, i, spp, sppl, sppr)) {
                distances[i] = 1;
            } else if (gNodeSNodeLinked(j, i, sppl, sppr)) {
                distances[i] = 2;
            } else {
                distances[i] = Integer.MAX_VALUE;
            }
        }
        for (int i = 0; i < gTree.getNodeCount(); i++) {
            if (distances[i] < Integer.MAX_VALUE) {
                while (!gTree.getNode(i).isRoot()) {
                    int anc = gTree.getNode(i).getParent().getNr();
                    distances[anc] = Math.min(distances[anc], distances[i] + 1);
                    i = anc;
                }
            }

        }
        fillInAnySubtreeDistances(gTree.getRoot(), distances);
    }



    private void fillInGTreeLogHeights(int j, double[] logHeights) {
        fillInAnyTreeLogHeights(gTrees.get(j), logHeights);
    }



    private void getGTreeFocusedScalerInternalBounds(double[] interval, int j,
                                                double[] logHeights, double[] weights) {
        getAnyTreeInternalFSBounds(interval, gTrees.get(j), logHeights, weights);
    }



    private void getGTreeFocusedScalerCompatibilityBounds(double[] interval, int j, OpFSinfoSMCNode fsNodeInfo,
                                                     double[] logHeights, double[] weights) {
        BitUnion sppL = fsNodeInfo.getLftUnion();
        BitUnion sppR = fsNodeInfo.getRgtUnion();
        double sLogHeight = fsNodeInfo.getLogHeight();
        double sWeight = fsNodeInfo.getWeight();
        ArrayList<Node> straddlers = unionArrays.getStraddlers(j, sppL, sppR);
        for (Node s : straddlers) {
            double weightDiff = sWeight - weights[s.getNr()];
            double loghDiff = logHeights[s.getNr()] - sLogHeight;
            if (weightDiff > 0) {
                interval[1] = Math.min(interval[1], loghDiff/weightDiff);
            }
            if (weightDiff < 0) {
                interval[0] = Math.max(interval[0], loghDiff/weightDiff);
            }
        }
    }



    private boolean gNodeSNodeCloselyLinked(int j, int i, BitUnion spp, BitUnion sppL, BitUnion sppR) {
        TreeInterface gTree = gTrees.get(j);
        Node iN = gTree.getNode(i);
        if (!iN.isLeaf()  &&  unionArrays.gNodeUnion(j, i).isContainedIn(spp)) {
            int lft = iN.getChild(0).getNr();
            int rgt = iN.getChild(1).getNr();
            boolean LwithinL = unionArrays.gNodeUnion(j, lft).isContainedIn(sppL);
            boolean RwithinL = unionArrays.gNodeUnion(j, rgt).isContainedIn(sppL);
            boolean LwithinR = unionArrays.gNodeUnion(j, lft).isContainedIn(sppR);
            boolean RwithinR = unionArrays.gNodeUnion(j, rgt).isContainedIn(sppL);
            return ((LwithinL & RwithinR) | (LwithinR & RwithinL));
        } else {
            return false;
        }
    }


    private boolean gNodeSNodeLinked(int j, int i, BitUnion sppL, BitUnion sppR) {
        TreeInterface gTree = gTrees.get(j);
        Node iN = gTree.getNode(i);
        if (!iN.isLeaf()) {
            boolean overlapsL = unionArrays.gNodeUnion(j, i).overlaps(sppL);
            boolean overlapsR = unionArrays.gNodeUnion(j, i).overlaps(sppR);
            boolean linked = overlapsL & overlapsR;
            int lft = iN.getChild(0).getNr();
            int rgt = iN.getChild(1).getNr();
            boolean LoverlapsL = unionArrays.gNodeUnion(j, lft).overlaps(sppL);
            boolean LoverlapsR = unionArrays.gNodeUnion(j, lft).overlaps(sppR);
            boolean Llinked = LoverlapsL & LoverlapsR;
            boolean RoverlapsL = unionArrays.gNodeUnion(j, rgt).overlaps(sppL);
            boolean RoverlapsR = unionArrays.gNodeUnion(j, rgt).overlaps(sppR);
            boolean Rlinked = RoverlapsL & RoverlapsR;
            return (linked & !Llinked & !Rlinked); // Findbugs doesn't like this but I do.
        } else {
            return false;
        }
    }



    /*****************************************************************************************/
    /**************************** Generic to SMC and gene trees ******************************/
    /*****************************************************************************************/

    private void fillInAnyTreeLogHeights(TreeInterface tree, double[] logHeights) {
        assert tree.getNodeCount() == logHeights.length;
        for (int i = 0; i < logHeights.length; i++) {
            if (!tree.getNode(i).isLeaf()) {
                logHeights[i] = Math.log(tree.getNode(i).getHeight());
            } else {
                logHeights[i] = Double.NEGATIVE_INFINITY;
            }
        }
    }



    private void fillInAnySubtreeDistances(Node node, int[] distances) {
        if (!node.isRoot()) {
            int anc = node.getParent().getNr();
            distances[node.getNr()] = Math.min(distances[node.getNr()], distances[anc] + 1);
        }
        if (!node.isLeaf()) {
            fillInAnySubtreeDistances(node.getChild(0), distances);
            fillInAnySubtreeDistances(node.getChild(1), distances);
        }
    }


    private double distanceToWeight(int distance, int rootDistance) {
        double alpha = 8.0;
        if (distance >= rootDistance) {
            return 0.0;
        } else {
            double x = (double)(rootDistance - distance) / (double)rootDistance;
            return (Math.pow(alpha, x) - 1.0) / (alpha - 1.0);
        }
    }


    // h() is logHeight, w() is weight. n is node, a is parent. x is amount to move.
    // This must be true after move:
    // h(a) + w(a)*x >= h(n) + w(n)*x
    // h(a) - h(n) >= (w(n) - w(a))x
    // If w(n) - w(a) > 0, get x <= (h(a) - h(n))/(w(n) - w(a))
    // If w(n) - w(a) < 0, get x >= (h(a) - h(n))/(w(n) - w(a))
    public void getAnyTreeInternalFSBounds(double[] interval, TreeInterface tree, double[] logHeights, double[] weights) {
        for (int n = 0; n < tree.getNodeCount(); n++) {
            if (!tree.getNode(n).isLeaf()  &&  !tree.getNode(n).isRoot()) {
                Node ancN = tree.getNode(n).getParent();
                double weightDiff = weights[n] - weights[ancN.getNr()];
                double loghDiff = logHeights[ancN.getNr()] - logHeights[n];
                if (weightDiff > 0) {
                    interval[1] = Math.min(interval[1], loghDiff/weightDiff);
                }
                if (weightDiff < 0) {
                    interval[0] = Math.max(interval[0], loghDiff/weightDiff);
                }
            }
        }
    }


    public double fsMoveDoScale(TreeInterface tree, double logSF, double [] weights) {
        double logHR = 0.0;
        for (int s = 0; s < tree.getNodeCount(); s++) {
            Node sN = tree.getNode(s);
            if (!sN.isLeaf() && !sN.isRoot() && weights[s] > 0.0) {
                double logSFW = logSF * weights[s];
                sN.setHeight(sN.getHeight() * Math.exp(logSFW));
                logHR += logSFW;

            }
        }
        return logHR;
    }

}

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

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


/**
 * Created by Graham Jones on 19/08/2014.
 */

@Description("A move which changes a node height in the SMC-tree, " +
             "and changes some node heights in the gene trees in a way which preserves all " +
             "topologies and usually maintains compatibility.")

// This is a multiple-trees operator. Is it a TreeOperator?
public class NodesNudge extends Operator {

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

    public enum GtreeNodeCriterion {
        MIXED_WITH_PURE_CHILDREN,
        MIXED_WITH_A_PURE_CHILD,
        MIXED
    }
    private final static double [] gtreeNodeCriterionWeights = {1,1,1};



    private Tree sTree;
    private List<Tree> gTrees;
    private long delay = 0;
    private int callCount = 0;
    private boolean sTreeTooSmall;

    private UnionArrays unionArrays;

    private int [] rejectCounts;
    private int [] acceptCounts;
    private double [] upTotals;
    private double [] downTotals;
    private int lastcriterion;
    private double lastheightchange;

    private final boolean debugFlag = Boolean.valueOf(System.getProperty("stacey.debug"));
    private int numberofdebugchecks = 0;
    private final static int maxnumberofdebugchecks = 100000;




    private static class OpNNinfoSMCNode {
        private final int nodeNr;
        private final double nodeHeight;
        private final double lftHeight;
        private final BitUnion nodeUnion;
        private final BitUnion lftUnion;
        private final double rgtHeight;
        private final BitUnion rgtUnion;
        private final double ancHeight;

        public OpNNinfoSMCNode(int nodeNr,
                               double nodeHeight, BitUnion nodeUnion,
                               double lftHeight, BitUnion lftUnion,
                               double rgtHeight, BitUnion rgtUnion,
                               double ancHeight) {
            this.nodeNr = nodeNr;
            this.nodeHeight = nodeHeight;
            this.nodeUnion = nodeUnion;
            this.lftHeight = lftHeight;
            this.lftUnion = lftUnion;
            this.rgtHeight = rgtHeight;
            this.rgtUnion = rgtUnion;
            this.ancHeight = ancHeight;
        }

        public int getNodeNr() { return nodeNr; }
        public double getNodeHeight() { return nodeHeight; }
        public BitUnion getNodeUnion() { return nodeUnion; }
        public double getLftHeight() { return lftHeight; }
        public BitUnion getLftUnion() { return lftUnion; }
        public double getRgtHeight() { return rgtHeight; }
        public BitUnion getRgtUnion() { return rgtUnion; }
        public double getAncHeight() { return ancHeight; }
    }



    /********************************************************************************************/


    @Override
    public void initAndValidate() {
        //super.initAndValidate(); // ??

        sTree = smcTreeInput.get();
        gTrees = geneTreesInput.get();
        sTreeTooSmall = (sTree.getLeafNodeCount() < 3);
        if (delayInput != null  &&  delayInput.get() != null) {
            delay = delayInput.get().longValue();
        }
        Bindings bindings = Bindings.initialise(sTree, gTrees);
        unionArrays = UnionArrays.initialise(sTree, gTrees, bindings);

        // debug and tuning stuff
        rejectCounts = new int[GtreeNodeCriterion.values().length];
        acceptCounts = new int[GtreeNodeCriterion.values().length];
        upTotals = new double[GtreeNodeCriterion.values().length];
        downTotals = new double[GtreeNodeCriterion.values().length];
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
            Checks.allTreesAndCompatibility(sTree, gTrees, "NodesNudge", "before move");
            numberofdebugchecks++;
        }

        unionArrays.update();
        double logHR = doNodesNudgeMove();     //  The business
        unionArrays.reset();

        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            Checks.allTreesAndCompatibility(sTree, gTrees, "NodesNudge", "after move");
            numberofdebugchecks++;
        }
        return logHR;
    }



    @Override
    public void accept() {
        super.accept();
        acceptCounts[lastcriterion]++;
        if (lastheightchange > 0) {
            upTotals[lastcriterion] += lastheightchange;
        } else {
            downTotals[lastcriterion] += lastheightchange;
        }


    }


    @Override
    public void reject(final int reason) {
        super.reject(reason);
        rejectCounts[lastcriterion]++;

    }


    @Override
    public void storeToFile(final PrintWriter out) {
        super.storeToFile(out);
        if (debugFlag) {
            out.print("{id:\"" + getID() + "\"");
            for (int ch = 0; ch < GtreeNodeCriterion.values().length;ch++) {
                out.print(" " + GtreeNodeCriterion.values()[ch].toString() +
                                ", accept:" + acceptCounts[ch] +
                                ", reject:" + rejectCounts[ch] +
                                ", downTotal:" + downTotals[ch] +
                                ", upTotals:" + upTotals[ch]
                );
            }
            out.print("}");
        }
    }


    /*************************************************************************************************/



    private double doNodesNudgeMove() {

        sTree.startEditing(this);
        // TODO This is a fix for beast.core.State$Trie memory
        for (int j = 0; j < gTrees.size(); j++) {
            gTrees.get(j).startEditing(this);
        }
        // what type of move
        assert GtreeNodeCriterion.values().length == gtreeNodeCriterionWeights.length;
        int choice = Randomizer.randomChoicePDF(gtreeNodeCriterionWeights);
        lastcriterion = choice;
        GtreeNodeCriterion criterion = GtreeNodeCriterion.values()[choice];
        return doNormalNodesNudgeMove(criterion);
    }



    // For doNodesNudgeMove(). All Nodes Nudge Moves are currently normal
    private double doNormalNodesNudgeMove(GtreeNodeCriterion criterion) {
        double logHR = 0.0;
        double [] interval;
        ArrayList< ArrayList<Node>> geneNodesToNudge = new ArrayList< ArrayList<Node>>(gTrees.size());

        // Choose the node in smcTree
        OpNNinfoSMCNode nani = chooseNodeForNudge();

        // Collect all nodes in all gene tree that are to be nudged.
        // TODO-threaded
        for (int j = 0; j < gTrees.size(); j++) {
            ArrayList<Node> pgtnodes = nodesToNudgeInGtree(j, nani, criterion);
            geneNodesToNudge.add(pgtnodes);
        }

        // get the constraints from gene trees
        interval = getMaxNodesNudgeIntervalForGeneTrees(geneNodesToNudge);

        // get the constraints from smcTree.
        interval[0] = Math.max(interval[0], nani.getLftHeight() - nani.getNodeHeight());
        interval[0] = Math.max(interval[0], nani.getRgtHeight() - nani.getNodeHeight());
        if (nani.getAncHeight() == Double.POSITIVE_INFINITY) {
            interval[1] = Math.min(interval[1], 0.1); // in case node is root and no other constraint
        } else {
            interval[1] = Math.min(interval[1], nani.getAncHeight() - nani.getNodeHeight());
        }
        assert interval[0] <= 0.0;
        assert interval[1] >= 0.0;

        // choose the amount to nudge by
        double dh = Randomizer.uniform(interval[0], interval[1]);
        logHR += doSTreeNodeNudge(nani.getNodeNr(), dh);
        lastheightchange = dh;

        logHR += doNodesNudgeOnGeneTrees(dh, geneNodesToNudge);
        return logHR;
    }



    // For doNormalNodesNudgeMove()
    private double[] getMaxNodesNudgeIntervalForGeneTrees(ArrayList< ArrayList<Node>> toNudge) {
        double [] interval = new double[2];
        interval[0] = -Double.MAX_VALUE;
        interval[1] = Double.MAX_VALUE;
        // TODO-threaded
        for (int pgt = 0; pgt < gTrees.size(); pgt++) {
            double [] pgtint = getMaxNodesNudgeIntervalWithinOneGeneTree(toNudge.get(pgt));
            interval[0] = Math.max(interval[0], pgtint[0]);
            interval[1] = Math.min(interval[1], pgtint[1]);
        }
        return interval;
    }



    // For getMaxNodesNudgeIntervalForGeneTrees()
    private double[] getMaxNodesNudgeIntervalWithinOneGeneTree(ArrayList<Node> nodes) {
        ArrayList<ArrayList<Node>> conncpts = findConnectedComponents(nodes);
        double[] interval = new double[2];
        interval[0] = -Double.MAX_VALUE;
        interval[1] = Double.MAX_VALUE;
        for (ArrayList<Node> cpt : conncpts) {
            Node root = cpt.get(0); // findConnectedComponentFrom() puts root first
            Node anc = root.getParent();
            if (anc != null) {
                interval[1] = Math.min(interval[1], anc.getHeight() - root.getHeight());
            }
            for (Node node : cpt) {
                double nhgt = node.getHeight();
                Node ch0 = node.getChild(0);
                Node ch1 = node.getChild(1);
                if (!cpt.contains(ch0)) {
                    interval[0] = Math.max(interval[0], ch0.getHeight() - nhgt);
                }
                if (!cpt.contains(ch1)) {
                    interval[0] = Math.max(interval[0], ch1.getHeight() - nhgt);
                }
            }
        }
        return interval;
    }



    // For getMaxNodesNudgeIntervalWithinOneGeneTree()
    private ArrayList<ArrayList<Node>> findConnectedComponents(ArrayList<Node> nodes) {
        ArrayList<ArrayList<Node>> conncpts = new ArrayList<ArrayList<Node>>();
        HashMap<Node, Integer> node2Index = new HashMap<Node, Integer>(nodes.size());
        int[] idxs = new int[nodes.size()];
        for (int i = 0; i < nodes.size(); i++) {
            node2Index.put(nodes.get(i), i);
        }
        for (int i = 0; i < nodes.size(); i++) {
            idxs[i] = -1;
        }
        int idx = 0;
        for (Node x : nodes) {
            // if not dealt with x yet, trace back to root of its connected component
            if (idxs[node2Index.get(x)] < 0) {
                while (x.getParent() != null && node2Index.containsKey(x.getParent())) {
                    x = x.getParent();
                }
                // find component from root and add it
                ArrayList<Node> cpt = findConnectedComponentFrom(x, node2Index);
                conncpts.add(cpt);
                // mark nodes in component as done (I'm only using idxs as flags)
                for (Node c : cpt) {
                    if (node2Index.get(c) == null) {
                        System.out.println("BUG in findConnectedComponents()");
                    }
                    assert node2Index.get(c) != null;
                    idxs[node2Index.get(c)] = idx;
                }
                idx++;
            }
        }
        return conncpts;
    }



    // For findConnectedComponents()
    private ArrayList<Node> findConnectedComponentFrom(Node x, HashMap<Node, Integer> node2Index) {
        ArrayList<Node> cpt = new ArrayList<Node>();
        cpt.add(x); // Note root added first
        // recursion from x within the set of nodes
        if (x.getChildCount() > 0) {
            Node ch0 = x.getChild(0);
            Node ch1 = x.getChild(1);
            if (ch0 != null && node2Index.containsKey(ch0)) {
                ArrayList<Node> cpt0 = findConnectedComponentFrom(ch0, node2Index);
                cpt.addAll(cpt0);
            }
            if (ch1 != null && node2Index.containsKey(ch1)) {
                ArrayList<Node> cpt1 = findConnectedComponentFrom(ch1, node2Index);
                cpt.addAll(cpt1);
            }
        }
        return cpt;
    }


    // iterates over gene trees to change heights
    private double doNodesNudgeOnGeneTrees(double dh, ArrayList< ArrayList<Node>> toNudge) {
        double logHR = 0.0;
        // TODO-threaded
        for (int j = 0; j < gTrees.size(); j++) {
            logHR += nudgeNodesInListBy(toNudge.get(j), dh);
        }
        return logHR;
    }


    // For doNodesNudgeOnGeneTrees(). This one actually edits the BEAST tree
    private double nudgeNodesInListBy(ArrayList<Node> nodes, double dh) {
        double logHR = 0.0;
        for (Node node : nodes) {
            node.setHeight(node.getHeight() + dh);
        }
        return logHR;
    }



    /**************************************************************************************/
    /*********************      SMC tree      ***************************************/
    /**************************************************************************************/

    public OpNNinfoSMCNode chooseNodeForNudge() {
        int s;
        do {
            s = Randomizer.nextInt(sTree.getNodeCount());
        } while (sTree.getNode(s).isLeaf()  ||  sTree.getNode(s).isRoot());
        Node sN = sTree.getNode(s);
        Node lftN = sN.getChild(0);
        Node rgtN = sN.getChild(1);
        Node ancN = sN.getParent();
        double sHgt = sN.getHeight();
        BitUnion sUnion = unionArrays.sNodeUnion(s);
        BitUnion lftUnion = unionArrays.sNodeUnion(lftN.getNr());
        BitUnion rgtUnion = unionArrays.sNodeUnion(rgtN.getNr());
        double lftHgt = lftN.getHeight();
        double rgtHgt = rgtN.getHeight();
        double anchgt = (ancN != null) ? ancN.getHeight() : Double.POSITIVE_INFINITY;
        return new OpNNinfoSMCNode(s, sHgt, sUnion, lftHgt, lftUnion, rgtHgt, rgtUnion, anchgt);
    }


    public double doSTreeNodeNudge(int s, double dh) {
        Node sN = sTree.getNode(s);
        sN.setHeight(sN.getHeight() + dh);
        return 0.0;
    }


    /**************************************************************************************/
    /*********************      gene trees      ***************************************/
    /**************************************************************************************/


    public ArrayList<Node> nodesToNudgeInGtree(int j, OpNNinfoSMCNode nani,
                                               NodesNudge.GtreeNodeCriterion criterion) {
        TreeInterface gTree = gTrees.get(j);
        BitUnion spp = nani.getNodeUnion();
        BitUnion sppl = nani.getLftUnion();
        BitUnion sppr = nani.getRgtUnion();
        ArrayList<Node> nodes = new ArrayList<Node>(0);
        for (int i = 0; i < gTree.getNodeCount(); i++) {
            Node iN = gTree.getNode(i);
            if (!iN.isLeaf()) {
                BitUnion gunion = unionArrays.gNodeUnion(j, i);
                if (gunion.isContainedIn(spp)) {
                    boolean mixed_with_pure_children = mixedWithPureChildren(j, i, sppl, sppr);
                    boolean mixed_with_a_pure_child = mixedWithAPureChild(j, i, sppl, sppr);
                    boolean mixed = mixedNode(j, i, sppl, sppr);
                    assert mixed || !mixed_with_pure_children;
                    assert mixed || !mixed_with_a_pure_child;
                    assert mixed_with_a_pure_child || !mixed_with_pure_children;

                    if (nodeWanted(criterion,
                            mixed_with_pure_children,
                            mixed_with_a_pure_child,
                            mixed)) {
                        nodes.add(iN);
                    }
                }
            }
        }
        return nodes;
    }




    // next 3 methods test various topological properties of nodes near smc tree node
    // which has tips in spp, and whose children have tips in sppl, sppr

    // this one is first meetings
    private boolean mixedWithPureChildren(int j, int i, BitUnion sppL, BitUnion sppR) {
        TreeInterface gTree = gTrees.get(j);
        int lft = gTree.getNode(i).getChild(0).getNr();
        int rgt = gTree.getNode(i).getChild(1).getNr();
        boolean lftInL = unionArrays.gNodeUnion(j, lft).isContainedIn(sppL);
        boolean rgtInL = unionArrays.gNodeUnion(j, rgt).isContainedIn(sppL);
        boolean lftInR = unionArrays.gNodeUnion(j, lft).isContainedIn(sppR);
        boolean rgtInR = unionArrays.gNodeUnion(j, rgt).isContainedIn(sppR);
        return ((lftInL && rgtInR) || (lftInR && rgtInL));
    }


    // this is a meeting where one child may be mixed but not both
    private boolean mixedWithAPureChild(int j, int i, BitUnion sppL, BitUnion sppR) {
        TreeInterface gTree = gTrees.get(j);
        BitUnion gunion = unionArrays.gNodeUnion(j, i);
        if (gunion.overlaps(sppL) && gunion.overlaps(sppR)) {
            int lft = gTree.getNode(i).getChild(0).getNr();
            int rgt = gTree.getNode(i).getChild(1).getNr();
            boolean lftInL = unionArrays.gNodeUnion(j, lft).isContainedIn(sppL);
            boolean rgtInL = unionArrays.gNodeUnion(j, rgt).isContainedIn(sppL);
            boolean lftInR = unionArrays.gNodeUnion(j, lft).isContainedIn(sppR);
            boolean rgtInR = unionArrays.gNodeUnion(j, rgt).isContainedIn(sppR);
            return (lftInL || rgtInR || lftInR || rgtInL);
        }
        return false;
    }


    // any meeting
    private boolean mixedNode(int j, int i, BitUnion sppL, BitUnion sppR) {
        BitUnion gunion = unionArrays.gNodeUnion(j, i);
        return (gunion.overlaps(sppL) && gunion.overlaps(sppR));
    }


    // various combinations of properties defined in last 5 methods
    private boolean nodeWanted(NodesNudge.GtreeNodeCriterion criterion,
                               boolean mixed_with_pure_children,
                               boolean mixed_with_a_pure_child,
                               boolean mixed) {
        boolean wanted;
        switch (criterion) {
            case MIXED_WITH_PURE_CHILDREN:
                wanted = mixed_with_pure_children;
                break;

            case MIXED_WITH_A_PURE_CHILD:
                wanted = mixed_with_a_pure_child;
                break;

            case MIXED:
                wanted = mixed;
                break;

            default:
                wanted = false;
                assert false;
        }
        return wanted;
    }



}




















/************************** Some experiments. See OperatorsGeneTree too **************************/

/*    public enum GtreeNodeCriterion {
        MIXED_WITH_PURE_CHILDREN_AND_PARENT_IN_GROUP,
        MIXED_WITH_PURE_CHILDREN,
        MIXED_WITH_A_PURE_CHILD,
        MIXED,
        PURE_MIXED_PARENT,
        PURE_MIXED_PARENT_OR_MIXED_WITH_PURE_CHILDREN,
        PURE_MIXED_PARENT_OR_MIXED_WITH_A_PURE_CHILD,
        PURE_MIXED_PARENT_OR_MIXED,
    }
    private final static double [] gtreeNodeCriterionWeights = {0,1,1,1,0,1,1,1};*/
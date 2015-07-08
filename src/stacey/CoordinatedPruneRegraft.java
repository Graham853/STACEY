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
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Created by Graham Jones on 27/10/2014.
 */

@Description("A move which does a subtree-prune-regraft move on a smcTree node," +
        "then carries a coordinated set of similar moves on the gene trees to maintain compatibility.")

// This is a multiple-trees operator. Is it a TreeOperator?
public class CoordinatedPruneRegraft extends Operator {

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<Tree> smcTreeInput =
            new Input<Tree>("smcTree",
                    "The species tree or minimal clusters tree", Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<List<Tree>> geneTreesInput =
            new Input<List<Tree>>("geneTree",
                    "All gene trees",
                    new ArrayList<Tree>());


    private TreeInterface sTree;
    private List<Tree> gTrees;
    private Bindings bindings;
    private UnionArrays unionArrays;


    private final boolean debugFlag = Boolean.valueOf(System.getProperty("stacey.debug"));
    private int numberofdebugchecks = 0;
    private final static int maxnumberofdebugchecks = 100000;
    private final int debugMaxBranches = 100;
    private int [] debugRejectCounts;
    private int [] debugAcceptCounts;
    private int debugLastNofBranches;



    // data for moving one subtree
    private class OpCPRinfoGTreeSpecification {
        private final Node source;
        private final double height;
        private final int choiceCount;
        private final Node destination;

        public OpCPRinfoGTreeSpecification(Node source, double height, int choiceCount, Node destination) {
            this.source = source;
            this.height = height;
            this.choiceCount = choiceCount;
            this.destination = destination;
        }


        public Node getSource() { return source; }
        public double getHeight() { return height; }
        public int getChoiceCount() { return choiceCount; }
        public Node getDestination() { return destination; }
    }


    private class OpCPRinfoSMCTreeSpecification {
        private final BitUnion sppS;
        private final BitUnion sppDs[];
        private final double heightsDancs[];

        public OpCPRinfoSMCTreeSpecification(BitUnion sppS,
                                             BitUnion[] sppDs,
                                             double[] heightsDancs) {
            this.sppS = sppS;
            this.sppDs = sppDs;
            this.heightsDancs = heightsDancs;
        }

        public BitUnion getSppS() { return sppS; }
        public BitUnion [] getSppDs() { return sppDs; }
        public double [] getHeightsDancs() { return heightsDancs; }
    }




    @Override
    public void initAndValidate() throws Exception {
        /*if (smcTreeInput.get().getLeafNodeCount() < 3) {
            throw new Exception("CoordinatedPruneRegraft cannot be used if there are less than 3 minimal clusters.");
        } TODO: this upsets Beauti */
        sTree = smcTreeInput.get();
        gTrees = geneTreesInput.get();
        bindings = Bindings.initialise(sTree, gTrees);
        unionArrays = UnionArrays.initialise(sTree, gTrees, bindings);
        if (debugFlag) {
            debugRejectCounts = new int[debugMaxBranches +1];
            debugAcceptCounts = new int[debugMaxBranches +1];
        }
    }



    @Override
    public double proposal() {

        if (smcTreeInput.get().getLeafNodeCount() < 3) {
            return Double.NEGATIVE_INFINITY;
        }

        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            Checks.allTreesAndCompatibility(sTree, gTrees, "CoordinatedPruneRegraft", "before move");
            numberofdebugchecks++;
        }

        unionArrays.update();
        double logHR = doCoordinatedPruneRegraft(false);// TODO: choose NNI and not-NNI at random(?)
        unionArrays.reset();

        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            Checks.allTreesAndCompatibility(sTree, gTrees, "CoordinatedPruneRegraft", "after move");
            numberofdebugchecks++;
        }
        return logHR;
    }


    @Override
    public void accept() {
        super.accept();
        if (debugFlag) {
            if (debugLastNofBranches < debugMaxBranches) {
                debugAcceptCounts[debugLastNofBranches]++;
            }
        }
    }


    @Override
    public void reject(final int reason) {
        super.reject(reason);
        if (debugFlag) {
            if (debugLastNofBranches < debugMaxBranches) {
                debugRejectCounts[debugLastNofBranches]++;
            }
        }
    }


    @Override
    public void storeToFile(final PrintWriter out) {
        super.storeToFile(out);
        if (debugFlag) {
            out.print("{id:\"" + getID() + "\"");
            for (int nb = 0; nb < debugMaxBranches; nb++) {
                out.print(" " + nb +
                                ", accept:" + debugAcceptCounts[nb] +
                                ", reject:" + debugRejectCounts[nb]
                );
            }
            out.print("}");
        }
    }



    // choose smc tree move s->d, remember s,x for reverse move, get information about
    // the move (smi) to pass to gene trees, do the smc tree move, do the gene tree
    // moves, get information about the reverse move (revsmi), calculate the
    // HR for reverse move.
    private double doCoordinatedPruneRegraft(boolean NNI) {

        // TODO This is a fix for beast.core.State$Trie memory
        for (int j = 0; j < gTrees.size(); j++) {
            gTrees.get(j).startEditing(null);
        }

        int d = -1;
        int s = -1;
        s = chooseSubtreeForSPR(NNI);
        Integer [] debugnb = new Integer[1];
        d = chooseDestinationForSPR(s, NNI, debugnb);
        debugLastNofBranches = debugnb[0];
        int x = siblingOf(sTree.getNode(s)).getNr();
        // Track where s and x end up by storing their unions (values not references!)
        BitUnion sUnion = new BitUnion(bindings.smcTreeTipCount());
        BitUnion xUnion = new BitUnion(bindings.smcTreeTipCount());
        sUnion.replaceWith(unionArrays.sNodeUnion(s));
        xUnion.replaceWith(unionArrays.sNodeUnion(x));
        // the move
        double logHR = 0.0;
        OpCPRinfoSMCTreeSpecification smi = makeSPRmoveInfo(s, d);
        doSPRmove(s, d);
        logHR += doAllCoordinatedSPRMoves(smi);
        // update the smcTree and gene trees, find s and x in new smcTree
        unionArrays.reset();
        unionArrays.update();
        int sNew = indexFromUnion(sUnion);
        int xNew = indexFromUnion(xUnion);
        // calculate the reverse HR
        OpCPRinfoSMCTreeSpecification revsmi = makeSPRmoveInfo(sNew, xNew);
        logHR -= HRForAllCoordinatedSPRMoves(revsmi);
        return logHR;
    }



    // Carries out a set of SPR moves on each gene tree based on the species tree SPR move
    // which has just been done. Also finds the Hastings ratio for the move.
    // This is log(number of choices) for the move.
    private double doAllCoordinatedSPRMoves(OpCPRinfoSMCTreeSpecification sppSPRInfo) {
        double logHR = 0;
        for (int j = 0; j < gTrees.size(); j++) {
            logHR += doCoordinatedSPRMoves(j, sppSPRInfo);
        }
        return  logHR;
    }


    // For calculation of Hastings ratio for SPR move. This finds log(number of choices)
    // for the reverse move.
    private double HRForAllCoordinatedSPRMoves(OpCPRinfoSMCTreeSpecification sppSPRInfo) {
        double logHR = 0;
        for (int j = 0; j < gTrees.size(); j++) {
            logHR += HRForCoordinatedSPRMoves(j, sppSPRInfo);
        }
        return  logHR;
    }


    /**********************************************************************************/
    /******************* Dealing with the SMC tree ************************************/
    /**********************************************************************************/

    /***************    CoordinatedPruneRegraft  ***************************************/
    /*
    1. Choose a node s (source). chooseSubtreeForSPR()
    2. Choose a node d (destination) so that d's branch contains time of s. chooseDestinationForSPR()
    3. Find info the gene trees need to know. makeSPRmoveInfo()
    4. Move the subtree with root s into branch d. doSPRmove()

         s                        .                         s
    x   /          d              .              x           \   d
     \ /          /               .               \           \ /
      p          /                .                \           p
       \        /                 .                 \         /
      --y      /                  .                --y       /
         \    z--e0               .                   \    z--
          \ /                     .                    \ /
           m = d2                 .                     m

    The gene trees need to know
    sppS the tips of s
    hM   the height of m, the mrca of s and d,
    sppD the tips of d
    sppEs the tips of the side branches between d and m, just e0 here
    heightsEancs the heights of the parents of these, just height of z here
    */

    // choose a random source node s, not root, not child of root.
    // If NNI, restrict to s with appropriate cousin.
    public int chooseSubtreeForSPR(boolean NNI) {
        assert sTree.getNodeCount() > 3;
        int s;
        if (NNI) {
            do {
                s = Randomizer.nextInt(sTree.getNodeCount());
            }
            while (rejectSubtreeForNNI(s));
        } else {
            do {
                s = Randomizer.nextInt(sTree.getNodeCount());
            }
            while (rejectSubtree(s));
        }
        return s;
    }



    private boolean rejectSubtree(int s) {
        return (sTree.getNode(s).isRoot()  ||  sTree.getNode(s).getParent().isRoot());
    }


    private boolean rejectSubtreeForNNI(int s) {
        Node sN = sTree.getNode(s);
        if (sN.isRoot()  ||  sN.getParent().isRoot()) {
            return true;
        }
        Node ancN = sN.getParent();
        Node sibancN = siblingOf(ancN);
        return  (sibancN.getHeight() >= ancN.getHeight());
    }



    // choose a destination for subtree s.  debugnb[1] is for debugging/tuning
    // Note: ancN is not root, so h < root height, so (dN.getHeight() <= h) implies dN not root
    // but maybe branch length of zero?
    private int chooseDestinationForSPR(int s, boolean NNI, Integer [] debugnb) {
        Node sN = sTree.getNode(s);
        if (NNI) {
            Node ancN = sN.getParent();
            debugnb[0] = 3;
            return siblingOf(ancN).getNr();
        } else {
            Node sibN = siblingOf(sN);
            Node ancN = sN.getParent();
            double h = ancN.getHeight();
            ArrayList<Integer> dests = new ArrayList<Integer>(1 + sTree.getNodeCount()/2);
            ArrayList<Integer> nbs = new ArrayList<Integer>(1 + sTree.getNodeCount()/2);
            for (int d = 0; d < sTree.getNodeCount(); d++) {
                Node dN = sTree.getNode(d);
                if (d != s  &&  d != ancN.getNr()  &&  d != sibN.getNr()) {
                    if (dN.getHeight() <= h  &&  !dN.isRoot()  &&  dN.getParent().getHeight() >= h) {
                        dests.add(d);
                        int nb = debugNumberOfBranchesBetween(s, d);
                        assert nb >= 3;
                        nbs.add(nb);
                    }
                }
            }
            int choice = Randomizer.nextInt(dests.size());
            int d = dests.get(choice);
            debugnb[0] = nbs.get(choice);
            return d;
        }
    }



    // Gather information that gene trees need to make moves corresponding to s->d
    private OpCPRinfoSMCTreeSpecification makeSPRmoveInfo(int s, int d) {
        BitUnion sppS = unionArrays.sNodeUnion(s);
        int m = mrcaOfPair(s, d);
        ArrayList<Integer> ds = new ArrayList<Integer>();
        while (d != m) {
            ds.add(d);
            d = sTree.getNode(d).getParent().getNr();
        }
        BitUnion [] sppDs = new BitUnion[ds.size()];
        double [] heightsDancs = new double[ds.size()];
        for (int i = 0; i < ds.size(); i++) {
            int di = ds.get(i);
            sppDs[i] = unionArrays.sNodeUnion(di);
            heightsDancs[i] = sTree.getNode(di).getParent().getHeight();
        }
        OpCPRinfoSMCTreeSpecification smi = new OpCPRinfoSMCTreeSpecification(sppS, sppDs, heightsDancs);
        return smi;
    }


    // this one actually edits the smc tree
    private void doSPRmove(int s, int d) {
        Node sr = sTree.getNode(s);
        Node dr = sTree.getNode(d);
        Node ar = sr.getParent();
        Node xr = ar.getChild(0) == sr ? ar.getChild(1) : ar.getChild(0);
        Node yr = ar.getParent();
        Node zr = dr.getParent();
        zr.removeChild(dr);
        ar.removeChild(xr);
        yr.removeChild(ar);
        zr.addChild(ar);
        ar.addChild(dr);
        yr.addChild(xr);
    }



    /****************************** low level routines for SMC tree *************************/

    // for finding indices of s and x after move so can calculate reverse HR
    private int indexFromUnion(BitUnion union) {
        return unionArrays.nodeIndexOfUnionInSubSTree(sTree.getRoot(), union);
    }




    private int mrcaOfPair(int x, int y) {
        while ( x != y) {
            Node xN = sTree.getNode(x);
            Node yN = sTree.getNode(y);
            if (xN.getHeight() < yN.getHeight()) {
                x = sTree.getNode(x).getParent().getNr();
            } else {
                y = sTree.getNode(y).getParent().getNr();
            }
        }
        return x;
    }


    private Node siblingOf(Node xN) {
        assert !xN.isRoot();
        Node ancN = xN.getParent();
        if (ancN.getChild(0) == xN) {
            return ancN.getChild(1);
        } else {
            return ancN.getChild(0);
        }
    }

    // for debugging
    private int debugNumberOfBranchesBetween(int x, int y) {
        int n = 0;
        while ( x != y) {
            Node xN = sTree.getNode(x);
            Node yN = sTree.getNode(y);
            if (xN.getHeight() < yN.getHeight()) {
                x = sTree.getNode(x).getParent().getNr();
            } else {
                y = sTree.getNode(y).getParent().getNr();
            }
            n++;
        }
        return n;

    }


    /**********************************************************************************/
    /******************* Dealing with gene trees **************************************/
    /**********************************************************************************/



    //  obtains a list of moves, sorts them, does them.
    private double doCoordinatedSPRMoves(int j, OpCPRinfoSMCTreeSpecification sppSPRInfo) {
        ArrayList<OpCPRinfoGTreeSpecification> gtsprss = makeListOfSPRSpecs(j, sppSPRInfo);
        Collections.sort(gtsprss, GTREESPRSPEC_ORDER);
        double logHR = 0;
        for (OpCPRinfoGTreeSpecification gtspr : gtsprss) {
            doSPRmove(gtspr.getSource(), gtspr.getDestination());
            logHR += Math.log(gtspr.getChoiceCount());
        }
        return  logHR;
    }


    //  obtains a list of moves, finds the HR, does not edit the tree. This is for the
    // reverse move to doCoordinatedSPRMoves().
    private double HRForCoordinatedSPRMoves(int j, OpCPRinfoSMCTreeSpecification sppSPRInfo) {
        ArrayList<OpCPRinfoGTreeSpecification> gtsprss = makeListOfSPRSpecs(j, sppSPRInfo);
        double logHR = 0;
        for (OpCPRinfoGTreeSpecification gtspr : gtsprss) {
            logHR += Math.log(gtspr.getChoiceCount());
        }
        return  logHR;
    }



    // makes a list of moves (GeneTreeSPRSpecifications) for one gene tree.
    // called from OperatorsGeneTree.
    private ArrayList<OpCPRinfoGTreeSpecification> makeListOfSPRSpecs(int j, OpCPRinfoSMCTreeSpecification sppSPRInfo) {
        TreeInterface gTree = gTrees.get(j);
        ArrayList<OpCPRinfoGTreeSpecification> gtsprss = new ArrayList<OpCPRinfoGTreeSpecification>(0);
        for (int a = 0; a < gTree.getNodeCount(); a++) {
            Node aN = gTree.getNode(a);
            double h = aN.getHeight();
            double [] heightsDancs = sppSPRInfo.getHeightsDancs();
            double hM = heightsDancs[heightsDancs.length-1];
            if (!aN.isLeaf()  &&  h <= hM) {
                int x = gTree.getNode(a).getChild(0).getNr();
                int y = gTree.getNode(a).getChild(1).getNr();
                boolean xinS = unionArrays.gNodeUnion(j, x).isContainedIn(sppSPRInfo.getSppS());
                boolean yinS = unionArrays.gNodeUnion(j, y).isContainedIn(sppSPRInfo.getSppS());
                if (xinS && !yinS) {
                    OpCPRinfoGTreeSpecification xspr =
                            chooseDestinationAndConstructMove(j, a, x, sppSPRInfo);
                    gtsprss.add(xspr);
                }
                if (yinS && !xinS) {
                    OpCPRinfoGTreeSpecification yspr =
                            chooseDestinationAndConstructMove(j, a, y, sppSPRInfo);
                    gtsprss.add(yspr);
                }
            }
        }
        return gtsprss;
    }




    // makes a list of possible destinations for one gene subtree
    // then chooses one and returns a GeneTreeSPRSpecification
    // (which is source, destination, height, number of choices).
    private OpCPRinfoGTreeSpecification chooseDestinationAndConstructMove(
                              int j, int anc, int s, OpCPRinfoSMCTreeSpecification sppSPRInfo) {
        TreeInterface gTree = gTrees.get(j);
        Node source = gTree.getNode(s);
        double height = gTree.getNode(anc).getHeight();

        double [] heightsDancs = sppSPRInfo.getHeightsDancs();
        int nDs = heightsDancs.length;
        assert height <= heightsDancs[nDs-1];
        int d = 0;
        while (height > heightsDancs[d]) {
            d++;
        }
        ArrayList<Node> possibleDestinations = new ArrayList<Node>();
        for (int x = 0; x < gTree.getNodeCount(); x++) {
            Node xN = gTree.getNode(x);
            if (xN.getHeight() <= height  &&  height <= xN.getParent().getHeight()  &&
                    unionArrays.gNodeUnion(j, x).isContainedIn(sppSPRInfo.getSppDs()[d])) {
                possibleDestinations.add(xN);
            }
        }
        int choiceCount = possibleDestinations.size();
        assert choiceCount > 0;
        Node destination = possibleDestinations.get(Randomizer.nextInt(choiceCount));
        OpCPRinfoGTreeSpecification spr = new OpCPRinfoGTreeSpecification(source, height, choiceCount, destination);
        return spr;
    }




    //  this one actually edits the BEAST gene tree
    private void doSPRmove(Node s, Node d) {
        Node a = s.getParent();
        Node x = a.getChild(0) == s ? a.getChild(1) : a.getChild(0);
        Node y = a.getParent();
        Node z = d.getParent();
        z.removeChild(d);
        a.removeChild(x);
        y.removeChild(a);
        z.addChild(a);
        a.addChild(d);
        y.addChild(x);
    }


    /****************************** low level routine for gene trees *************************/

    // used for sorting the SPR moves from most ancient to most recent.
    private static final Comparator<OpCPRinfoGTreeSpecification> GTREESPRSPEC_ORDER =
            new Comparator<OpCPRinfoGTreeSpecification>() {
                public int compare(OpCPRinfoGTreeSpecification a, OpCPRinfoGTreeSpecification b) {
                    if (a.getHeight() != b.getHeight()) {
                        if (a.getHeight() < b.getHeight()) {
                            return 1;
                        } else {
                            return -1;
                        }
                    } else {
                        return 0;
                    }
                }
            };

}

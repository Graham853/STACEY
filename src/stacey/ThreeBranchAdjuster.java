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
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import stacey.debugtune.Checks;
import stacey.util.Bindings;
import stacey.util.BitUnion;
import stacey.util.UnionArrays;
import stacey.BirthDeathCollapseModel;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import static stacey.BirthDeathCollapseModel.belowCollapseHeight;

/**
 * Created by Graham Jones on 09/08/2015.
 *
 *
 */

@Description("An implementation of the rubber-band operator from Yang and Rannala (2003)." +
             "It changes the height of a SMC-tree node, and all gene tree nodes in the three adjacent branches.")
public class ThreeBranchAdjuster extends Operator {


    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<Tree> smcTreeInput =
            new Input<Tree>("smcTree",
                    "The species tree or minimal clusters tree", Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<List<Tree>> geneTreesInput =
            new Input<List<Tree>>("geneTree",
                    "All gene trees",
                    new ArrayList<Tree>());

    public Input<RealParameter> popSFInput =
            new Input<>("popSF",
                    "The population scaling factor for the STACEY coalescent", Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<Double> tuningInput =
            new Input<>("tuning",
                    "A tuning parameter. The default is 1.0, and allowed values are positive numbers. " +
                    "Larger values make bigger jumps and so decrease the acceptance ratio. " +
                            "Experiments in the range [0.1,10.0] seem sensible.");

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<Long> delayInput =
            new Input<>("delay",
                    "Number of times the operator is disabled.");

    private Tree sTree;
    private List<Tree> gTrees;
    private int callCount = 0;
    private boolean sTreeTooSmall;

    private long delay = 0;
    private double tuning = 1.0;


    private UnionArrays unionArrays;

    private final boolean debugFlag = Boolean.valueOf(System.getProperty("stacey.debug"));
    private int numberofdebugchecks = 0;
    private final static int maxnumberofdebugchecks = 100000;




    @Override
    public void initAndValidate() {
        //super.initAndValidate(); // ??

        sTree = smcTreeInput.get();
        gTrees = geneTreesInput.get();
        sTreeTooSmall = (sTree.getLeafNodeCount() < 3);
        if (delayInput.get() != null  &&  delayInput.get() != null) {
            delay = delayInput.get().longValue();
        }
        Bindings bindings = Bindings.initialise(sTree, gTrees);
        unionArrays = UnionArrays.initialise(sTree, gTrees, bindings);
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
            Checks.allTreesAndCompatibility(sTree, gTrees, "ThreeBranchAdjuster", "before move");
            numberofdebugchecks++;
        }

        unionArrays.update();
        double logHR = doThreeBranchAdjust();     //  The business
        unionArrays.reset();

        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            Checks.allTreesAndCompatibility(sTree, gTrees, "ThreeBranchAdjuster", "after move");
            numberofdebugchecks++;
        }
        return logHR;
    }



    private double doThreeBranchAdjust() {
        double logHR = 0.0;

        // TODO This is a fix for beast.core.State$Trie memory
        for (int j = 0; j < gTrees.size(); j++) {
            gTrees.get(j).startEditing(null);
        }
        // choose sTree node
        int sNodeNr;
        do {
            sNodeNr = Randomizer.nextInt(sTree.getNodeCount());
        } while (sTree.getNode(sNodeNr).isLeaf());
        Node sNode = sTree.getNode(sNodeNr);
        boolean isRoot = sTree.getNode(sNodeNr).isRoot();
        // Info about three branches meeting sNode
        double sHeight = sNode.getHeight();
        double sAncHeight;
        if (isRoot) {
            double maxGtreeHeight = Double.NEGATIVE_INFINITY;
            for (int j = 0; j <gTrees.size(); j++) {
                Tree gTree = gTrees.get(j);
                maxGtreeHeight = Math.max(maxGtreeHeight, gTree.getRoot().getHeight());
            }
            sAncHeight = maxGtreeHeight;
        } else {
            sAncHeight = sNode.getParent().getHeight();
        }

        double sLftHeight = sNode.getChild(0).getHeight();
        double sRgtHeight = sNode.getChild(1).getHeight();
        BitUnion sUnion = unionArrays.sNodeUnion(sNodeNr);
        BitUnion sLftUnion = unionArrays.sNodeUnion(sNode.getChild(0).getNr());
        BitUnion sRgtUnion = unionArrays.sNodeUnion(sNode.getChild(1).getNr());

        // collect all gTree nodes in these branches
        ArrayList<Node> gNodes = new ArrayList<>();
        ArrayList<Node> gLftNodes = new ArrayList<>();
        ArrayList<Node> gRgtNodes = new ArrayList<>();
        // TODO-threaded -  the ArrayLists look tricky
        for (int j = 0; j <gTrees.size(); j++) {
            Tree gTree = gTrees.get(j);
            for (int i = gTree.getLeafNodeCount(); i < gTree.getNodeCount(); i++) {
                double yHeight = gTree.getNode(i).getHeight();
                BitUnion yUnion = unionArrays.gNodeUnion(j, i);
                int timesAdded = 0; // debug check
                if (yUnion.isContainedIn(sUnion)  &&  yHeight < sAncHeight  &&  yHeight >= sHeight) {
                    gNodes.add(gTree.getNode(i));
                    timesAdded ++;
                }
                if (yUnion.isContainedIn(sLftUnion)  &&  yHeight < sHeight  &&  yHeight > sLftHeight) {
                    gLftNodes.add(gTree.getNode(i));
                    timesAdded ++;
                }
                if (yUnion.isContainedIn(sRgtUnion)  &&  yHeight < sHeight  &&  yHeight > sRgtHeight) {
                    gRgtNodes.add(gTree.getNode(i));
                    timesAdded ++;
                }
                assert timesAdded <= 1;
            }
        }
        double minSH = Math.max(sLftHeight, sRgtHeight);
        double maxSH = sAncHeight;
        double popSF = popSFInput.get().getValue();
        if (tuningInput != null  &&  tuningInput.get() != null) {
            tuning =  tuningInput.get().doubleValue();
        }
        double hgtRange = Math.min(tuning * 40.0 * popSF / gTrees.size(), 0.5*(maxSH - minSH));
        double newHeight = sHeight + Randomizer.uniform(-hgtRange, hgtRange);
        if (newHeight < minSH) {
            newHeight = 2 * minSH - newHeight;
        } else if (newHeight > maxSH) {
            newHeight = 2 * maxSH - newHeight;
        }
        sNode.setHeight(newHeight);
        logHR += setNewGNodeHeightsFromSHeight(gNodes,    sHeight, newHeight, sAncHeight);
        logHR += setNewGNodeHeightsFromSHeight(gLftNodes, sHeight, newHeight, sLftHeight);
        logHR += setNewGNodeHeightsFromSHeight(gRgtNodes, sHeight, newHeight, sRgtHeight);

        /* non-uniform version, Jan 2016. Seems no better, perhaps worse, KISS applies. */
        /*
        double eps = BirthDeathCollapseModel.collapseHeight();
        double gap = eps - minSH;
        if (Double.compare(collapseWeightInput.get().getValue(), 0.0) == 0  ||
                gap <= 1e-8  ||
                maxSH - minSH <= 3.0 * gap) {
            // normal case
            double newHeight = sHeight + Randomizer.uniform(-hgtRange, hgtRange);
            if (newHeight < minSH) {
                newHeight = 2 * minSH - newHeight;
            } else if (newHeight > maxSH) {
                newHeight = 2 * maxSH - newHeight;
            }
            sNode.setHeight(newHeight);
            logHR += setNewGNodeHeightsFromSHeight(gNodes,    sHeight, newHeight, sAncHeight);
            logHR += setNewGNodeHeightsFromSHeight(gLftNodes, sHeight, newHeight, sLftHeight);
            logHR += setNewGNodeHeightsFromSHeight(gRgtNodes, sHeight, newHeight, sRgtHeight);
        } else {
            // doing species delimitation with oldest child 'sensibly' below collapse height
            // and enough room above (maxSH big enough) that the non-uniform sampling is sensible
            // The PDF has support [minSH, maxSH] and in this range is
            //  f(x) = [log(maxSH - minSH + delta) - log(delta)]^-1  /  (x - minSH + delta)
            //  The CDF is
            //  F(x) = [log(maxSH - minSH + delta) - log(delta)]^-1  *  [log(x - minSH + delta) - log(delta)]
            // The inverse CDF is
            // G(y) = minSH - delta + exp{  y * [log(maxSH - minSH + delta] - log(delta)) + log(delta)  }

            assert gap > 1e-8;
            double delta = gap*gap / (maxSH - minSH - 2*gap);
            assert delta > 0.0;
            double newHeight = newHeightSample(minSH, maxSH, delta);
            double oldDensity = newHeightPDF(sHeight, minSH, maxSH, delta);
            double newDensity = newHeightPDF(newHeight, minSH, maxSH, delta);
            logHR = Math.log(oldDensity/newDensity);
            sNode.setHeight(newHeight);
            logHR += setNewGNodeHeightsFromSHeight(gNodes,    sHeight, newHeight, sAncHeight);
            logHR += setNewGNodeHeightsFromSHeight(gLftNodes, sHeight, newHeight, sLftHeight);
            logHR += setNewGNodeHeightsFromSHeight(gRgtNodes, sHeight, newHeight, sRgtHeight);
        }

        Two functions used by this idea:
            private double newHeightPDF(double h, double minSH, double maxSH, double delta) {
                double d = 1.0 / (Math.log(maxSH - minSH + delta) - Math.log(delta));
                d /=  (h - minSH + delta);
                return d;
            }


            private double newHeightSample(double minSH, double maxSH, double delta) {
                double y = Randomizer.nextDouble();

                double q = minSH - delta +
                        Math.exp(  y * (Math.log(maxSH - minSH + delta) - Math.log(delta))+ Math.log(delta)  );
                assert q >= minSH - 1e-12;
                assert q <= maxSH + 1e-12;
                q = Math.max(q, minSH);
                q = Math.min(q, maxSH);
                return q;
            }


        */
        return logHR;
    }







    private double setNewGNodeHeightsFromSHeight(ArrayList<Node> nodes, double oldSH, double newSH, double limit) {
        double logHR = 0.0;
        // TODO-threaded. Not good, too fine-grained
        for (int i = 0; i < nodes.size(); i++) {
            Node node = nodes.get(i);
            double oldGH = node.getHeight();
            double scale = (newSH - limit) / (oldSH - limit);
            double newGH = limit + (oldGH - limit) * scale;
            node.setHeight(newGH);
            logHR += Math.log(scale);
        }
        return logHR;
    }



    private void setNewGNodeHeightsFromWeights(ArrayList<Node> nodes, ArrayList<Double> weights, double eta) {
        for (int i = 0; i < nodes.size(); i++) {
            Node node = nodes.get(i);
            double wt = weights.get(i+1);
            node.setHeight(node.getHeight() + wt * eta);
        }
    }


    private void updateBoundsFromArrays(double[] interval, double firstH, ArrayList<Node> nodes, double lastH, ArrayList<Double> weights) {
        int n = nodes.size();
        if (n > 0) {
            updateBounds(interval, firstH, nodes.get(0).getHeight(), weights.get(0), weights.get(1)) ;
            for (int i = 0; i+1 < n; i++) {
                updateBounds(interval, nodes.get(i).getHeight(), nodes.get(i+1).getHeight(), weights.get(i+1), weights.get(i+2));
            }
            updateBounds(interval, nodes.get(n-1).getHeight(), lastH, weights.get(n), weights.get(n+1));
        } else {
            updateBounds(interval, firstH, lastH, weights.get(0), weights.get(1));
        }
    }


    // h1 <= h2 are heights, w1, w2 are corresponding weights. This ensures that
    // h2 + w2*x >= h1 + w1*x for any x in (interval[0], interval[1])
    // (w2 - w1)*x >= h1 - h2
    // If w2 - w1 > 0, need x >= (h1 - h2)/(w2 - w1)
    // If w2 - w1 < 0, need x <= (h1 - h2)/(w2 - w1)
    private void updateBounds(double [] interval, double h1, double h2, double w1, double w2) {
        double wDiff = w2 - w1;
        double hDiff = h1 - h2;
        if (wDiff > 0) {
            interval[0] = Math.max(interval[0], hDiff/wDiff);
        }
        if (wDiff < 0) {
            interval[1] = Math.min(interval[1], hDiff/wDiff);
        }
    }



    // sorts the coalescences on height
    private static final Comparator<Node> COALESCENCE_ORDER = new Comparator<Node>() {
        public int compare(Node a, Node b) {
            return Double.compare(a.getHeight(), b.getHeight());
        }
    };

}

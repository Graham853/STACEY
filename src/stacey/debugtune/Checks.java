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

package stacey.debugtune;


import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import stacey.GtreeAndCoalFactor;
import stacey.util.*;

import java.util.ArrayList;
import java.util.List;


public class Checks {


    public static void allTreesAndCompatibility(TreeInterface sTree, List<Tree> gTrees,
                                                String opName, String when) {
        if (!Checks.treeIsOK(sTree)) {
            System.err.println("BUG found in " + opName + ". Bad sTree " + when);
            System.err.println(Misc.allTreesAsText(sTree, gTrees));
            throw new RuntimeException("Fatal STACEY error.");
        }
        for (TreeInterface gTree : gTrees) {
            if (!Checks.treeIsOK(gTree)) {
                System.err.println("BUG found in " + opName + ". Bad gTree " + when);
                System.err.println(Misc.allTreesAsText(sTree, gTrees));
                throw new RuntimeException("Fatal STACEY error.");
            }
        }
        if (!Checks.compatible(sTree, gTrees)) {
            System.err.println("BUG found in " + opName + ". Incompatible trees " + when);
            System.err.println(Misc.allTreesAsText(sTree, gTrees));
            throw new RuntimeException("Fatal STACEY error.");
        }
    }


    public static boolean allTreesAreOKAndCompatible(TreeInterface sTree, List<Tree> gTrees) {
        if (!Checks.treeIsOK(sTree)) {
            return false;
        }
        for (TreeInterface gTree : gTrees) {
            if (!Checks.treeIsOK(gTree)) {
                return false;
            }
        }
        if (!Checks.compatible(sTree, gTrees)) {
            return false;
        }
        return true;
    }



    public static boolean treeIsOK(TreeInterface tree) {
        Node nodes [] = tree.getNodesAsArray();
        for (int n = 0; n <nodes.length; n++) {
            if (!Double.isFinite(nodes[n].getHeight())) {
                return false;
            }
            if (nodes[n].getHeight() < 0) {
                return false;
            }
        }
        for (int n = 0; n <nodes.length; n++) {
            if (!nodes[n].isLeaf()) {
                if (nodes[n].getChildCount() != 2) {
                    return false;
                }
                double nodeHeight = nodes[n].getHeight();
                if (nodes[n].getChild(0).getHeight() > nodeHeight) {
                    return false;
                }
                if (nodes[n].getChild(1).getHeight() > nodeHeight) {
                    return false;
                }
            }
        }
        return true;
    }


    public static boolean compatible(TreeInterface sTree, List<Tree> gTrees) {
        Bindings bindings = Bindings.initialise(sTree, gTrees);
        UnionArrays unionArrays = UnionArrays.initialise(sTree, gTrees, bindings);
        unionArrays.update();
        if (!allUnionsLookOK(sTree, gTrees, unionArrays)) {
            return false;
        }
        for (int j = 0; j < gTrees.size(); j++) {
            TreeInterface gTree = gTrees.get(j);
            for (int i = 0; i < gTree.getNodeCount(); i++) {
                BitUnion gUnion = unionArrays.gNodeUnion(j, i);
                double gHeight = gTree.getNode(i).getHeight();
                int n = unionArrays.nodeIndexOfUnionInSubSTree(sTree.getRoot(), gUnion);
                if (sTree.getNode(n).getHeight() > gHeight) {
                    return false;
                }
            }

        }
        unionArrays.reset();
        return true;
    }


    private static boolean allUnionsLookOK(TreeInterface sTree, List<Tree> gTrees, UnionArrays unionArrays) {
        BitUnion [] sUnions = new BitUnion[sTree.getNodeCount()];
        for (int n = 0; n <sTree.getNodeCount(); n++) {
            sUnions[n] = unionArrays.sNodeUnion(n);
        }
        if (!unionsLookOK(sTree, sUnions)) {
            return false;
        }
        for (int j = 0; j < gTrees.size(); j++) {
            TreeInterface gTree = gTrees.get(j);
            BitUnion [] gUnions = new BitUnion[gTree.getNodeCount()];
            for (int i = 0; i < gTree.getNodeCount(); i++) {
                gUnions[i] = unionArrays.gNodeUnion(j, i);
            }
            if (!unionsLookOK(gTree, gUnions)) {
                return false;
            }
        }
        return true;
    }


    private static boolean unionsLookOK(TreeInterface tree, BitUnion [] unions) {
        for (int n = 0; n<tree.getNodeCount(); n++) {
            Node node = tree.getNode(n);
            if (node.isLeaf()) {
                BitUnion hereU = unions[node.getNr()];
                if (hereU.debugNumberBitsSets() != 1) {
                    return false;
                }
            } else {
                BitUnion hereU = unions[node.getNr()];
                BitUnion lftU  = unions[node.getChild(0).getNr()];
                BitUnion rgtU  = unions[node.getChild(1).getNr()];
                int size = hereU.debugSize();
                for (int b = 0; b < size; b++) {
                    boolean here = hereU.debugBitIsSet(b);
                    boolean lft = lftU.debugBitIsSet(b);
                    boolean rgt = rgtU.debugBitIsSet(b);
                    if ((lft | rgt) != here) {
                        return false;
                    }
                }
            }
            if (node.isRoot()) {
                BitUnion hereU = unions[node.getNr()];
                if (hereU.debugNumberBitsSets() != hereU.debugSize()) {
                    return false;
                }
            }
        }
        return true;
    }

}

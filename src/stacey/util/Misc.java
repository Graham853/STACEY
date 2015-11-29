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

package stacey.util;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;


import java.util.*;



@Description("Miscellaneous methods for STACEY, mainly for debug output.")
public class Misc {

    public static String allTreesAsText(TreeInterface sTree, List<Tree> gTrees) {
        StringBuffer buffer = new StringBuffer();
        buffer.append("***************************   All trees   ********************************\n");
        Bindings bindings = Bindings.initialise(sTree, gTrees);
        UnionArrays unionArrays = UnionArrays.initialise(sTree, gTrees, bindings);
        unionArrays.update();

        BitUnion [] sUnions = new BitUnion[sTree.getNodeCount()];
        for (int n = 0; n <sTree.getNodeCount(); n++) {
            sUnions[n] = unionArrays.sNodeUnion(n);
        }
        buffer.append(beastSubtreeAsTextWithHeader(sTree.getRoot(), sUnions));
        for (int j = 0; j < gTrees.size(); j++) {
            TreeInterface gTree = gTrees.get(j);
            BitUnion [] gUnions = new BitUnion[gTree.getNodeCount()];
            for (int i = 0; i < gTree.getNodeCount(); i++) {
                gUnions[i] = unionArrays.gNodeUnion(j, i);
            }
            buffer.append("\nGene tree" + j + "\n");
            buffer.append(beastSubtreeAsTextWithHeader(gTree.getRoot(), gUnions));
        }
    return buffer.toString();
    }



    public static String beastSubtreeAsNewick(TreeInterface tree, Node node) {
        ArrayList<String> nodeLabels = new ArrayList<>(tree.getNodeCount());
        for (int n = 0; n < tree.getNodeCount(); n++) {
            nodeLabels.add(tree.getNode(n).getID() == null ? "" + n : tree.getNode(n).getID());
        }
        return node.toString(nodeLabels);
    }


    public static String beastSubtreeAsTextWithHeader(Node node, BitUnion [] unions) {
        String header = "topology                                 idx height\n";
        String s = "";
        Stack<Integer> x = new Stack<Integer>();
        return header + beastSubtreeAsText(node, s, x, 0, "", unions);
    }


    private static String beastSubtreeAsText(Node node, String s, Stack<Integer> x, int depth, String b, BitUnion [] unions) {
        Integer[] y = x.toArray(new Integer[x.size()]);
        StringBuffer indent = new StringBuffer();
        for (int i = 0; i < depth; i++) {
            indent.append("  ");
        }
        for (int i = 0; i < y.length; i++) {
            indent.replace(2 * y[i], 2 * y[i] + 1, "|");
        }
        if (b.length() > 0) {
            indent.replace(indent.length() - b.length(), indent.length(), b);
        }
        s += indent;
        s += beastNodeAsText(node, indent.length(), unions);
        s += System.getProperty("line.separator");
        String subs = "";
        if (!node.isLeaf()) {
            x.push(depth);
            subs += beastSubtreeAsText(node.getChild(0), "", x, depth + 1, "-", unions);
            x.pop();
            subs += beastSubtreeAsText(node.getChild(1), "", x, depth + 1, "`-", unions);
        }
        return s + subs;
    }


    private static String beastNodeAsText(Node node, int indentlen, BitUnion [] unions) {
        StringBuilder s = new StringBuilder();
        Formatter formatter = new Formatter(s, Locale.US);
        if (node.isLeaf()) {
            formatter.format("%s ", node.getID());
        } else {
            formatter.format("%s ", "+");
        }
        while (s.length() < 40 - indentlen) {
            formatter.format("%s", " ");
        }
        int n = node.getNr();
        formatter.format("%3d ", n);
        formatter.format("%s ", nonnegIn8Chars(node.getHeight()));
        formatter.format("%30s", unions[node.getNr()].asText());
        return s.toString();
    }


    public static String nonnegIn8Chars(double x) {
        StringBuilder s = new StringBuilder();
        Formatter formatter = new Formatter(s, Locale.US);
        if (x < 0) {
            formatter.format("%8s", "NA");
        } else if (x == 0.0) {
            formatter.format("%8s", "zero");
        } else if (x < 1e-3) {
            formatter.format("%8.2e", x);
        } else if (x < 9.999) {
            formatter.format("%8.5f", x);
        } else if (x < 99.99) {
            formatter.format("%8.4f", x);
        } else if (x < 999.9) {
            formatter.format("%8.3f", x);
        } else if (x < 9999) {
            formatter.format("%8.2f", x);
        } else {
            formatter.format("%8.0f", x);
        }
        return s.toString();
    }

}
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

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Tree;


@Description("A gene tree and a coalescence factor (ploidy)")
public class GtreeAndCoalFactor extends CalculationNode {

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<Tree> gtree =
            new Input<Tree>("tree",
                    "A gene tree");

    @SuppressWarnings("CanBeFinal")
    public Input<Double> coalFactor =
            new Input<Double>("Ploidy",
                    "Coalescence factor, often the same as ploidy. Default is 2",
                    2.0);


    @Override
    public void initAndValidate() {


    }

    Tree getTree() {
        return gtree.get();
    }

    public double getCoalFactor() {
        return coalFactor.get();
    }
}

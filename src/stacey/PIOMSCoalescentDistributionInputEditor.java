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
// copied SpeciesTreePriorInputEditor then changed

package stacey;

import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.BEASTObject;
import beast.core.Input;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;



public class PIOMSCoalescentDistributionInputEditor extends InputEditor.Base  {

    public PIOMSCoalescentDistributionInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return PIOMSCoalescentDistribution.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
        super.init(input, plugin, itemNr, bExpandOption, bAddButtons);
    }

    protected void addComboBox(JComponent box, Input<?> input, BEASTObject plugin) {
        m_bAddButtons = true;
        String label = "Species Tree Population Size";
        addInputLabel(label, label);
        m_bAddButtons = false;
        add(Box.createHorizontalGlue());
    }

}

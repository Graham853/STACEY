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
// Copied from GeneTreeForSpeciesTreeDistributionInputEditor then edited

package stacey;

import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.BEASTObject;
import beast.core.Input;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;



public class GtreeAndCoalFactorInputEditor  extends InputEditor.Base {
    private static final long serialVersionUID = 1L;

    public GtreeAndCoalFactorInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return GtreeAndCoalFactor.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, InputEditor.ExpandOption bExpandOption, boolean bAddButtons) {
        m_bAddButtons = bAddButtons;
        m_input = input;
        m_plugin = plugin;
        this.itemNr= itemNr;
        String sID = plugin.getID();
        if (sID.contains(".t:")) {
            sID = sID.substring(sID.indexOf(".t:") + 3);
        }
        add(new JLabel("Gene Tree " + sID));
        add(Box.createGlue());
    }

    private static final int OTHER = 3;
    private final String [] sValues = new String[]{"autosomal_nuclear", "X", "Y or mitochondrial", "other"};
    private final Double [] fValues = new Double[]{2.0, 1.5, 0.5, -1.0};
    private JComboBox m_selectPluginBox;

    public InputEditor createPloidyEditor() {
        InputEditor editor = new InputEditor.Base(doc) {
            @Override
            public Class<?> type() {
                // TODO Auto-generated method stub
                return null;
            }

            @Override
            public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
                m_plugin = plugin;
                m_input = input;
                m_bAddButtons = bAddButtons;
                this.itemNr = itemNr;
                addInputLabel();

                // GRJ added <String> to avoid warning
                m_selectPluginBox = new JComboBox<String>(sValues);
                setSelection();
                String sSelectString = input.get().toString();
                m_selectPluginBox.setSelectedItem(sSelectString);

                m_selectPluginBox.addActionListener(new ActionListener() {
                    // implements ActionListener
                    public void actionPerformed(ActionEvent e) {
                        int i = m_selectPluginBox.getSelectedIndex();
                        if (i == OTHER) {
                            setSelection();
                            return;
                        }
                        try {
                            setValue(fValues[i]);
                            //lm_input.setValue(sSelected, m_plugin);
                        } catch (Exception e1) {
                            e1.printStackTrace();
                        }
                    }
                });
                m_selectPluginBox.setToolTipText(input.getHTMLTipText());
                add(m_selectPluginBox);
                add(Box.createGlue());
            }

            private void setSelection() {
                Double value = (Double) m_input.get();
                m_selectPluginBox.setSelectedIndex(OTHER);
                for (int i = 0; i < fValues.length; i++) {
                    if (value.equals(fValues[i])) {
                        m_selectPluginBox.setSelectedIndex(i);
                    }
                }
            }

        };
        editor.init(((GtreeAndCoalFactor)m_plugin).coalFactor,
                m_plugin, -1, InputEditor.ExpandOption.FALSE, true);
        return editor;
    }




}

package BDAStar;

import kstar.InteractionGraph;
import kstar.Molecule;
import kstar.PairwiseEnergyMatrix;

public class EnergyFunction 
{

    private PairwiseEnergyMatrix eMatrix;
    private InteractionGraph G;
    
    public EnergyFunction (PairwiseEnergyMatrix m, InteractionGraph g)
    {
    	eMatrix = m;
    	G = g;
    }

    public double score(ProteinConformation c)
    {
        return computeEforState(c);
    }

    private double computeEforState(ProteinConformation c)
    {

        int pi=0,ai=0,ri=0,pj=0,aj=0,rj=0;
        float energy = eMatrix.getShellShellE();
        /*
        en[0] =eMatrix.getShellShellE(); //Add shell shell energy
        int numPos = M.size() + lambda.size();
         */
        for(Position p: c.getPositions())
        {

            ProteinChoice choice = (ProteinChoice)c.getChoiceAt(p);
            ProteinPosition position = (ProteinPosition) p;
            pi=position.designIndex;
            ai=choice.aminoAcid;
            ri=choice.rotamer;

            /*
                en[0]+=eMatrix[pi][ai][ri][pi][0][1]; //add the self energy of the rotamer of the lambda residue
                en[0]+=eMatrix[pi][ai][ri][pi][0][0];
             */
            energy +=eMatrix.getShellRotE(pi, ai, ri);

            for(Position p2 : c.getPositions())
            {
                ProteinChoice choice2 = (ProteinChoice) c.getChoiceAt(p2);
                ProteinPosition position2 = (ProteinPosition) p2;
                
                pj=position2.designIndex;
                aj=choice2.aminoAcid;
                rj=choice2.rotamer;

                if(G.edgeExists(pi,pj))
                {  // if edge exists between residues in the interaction graph
                    energy +=eMatrix.getPairwiseE(pi, ai, ri, pj, aj, rj);
                }
            }
        }
        return energy;
    } 
}

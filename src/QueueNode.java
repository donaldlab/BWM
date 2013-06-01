///////////////////////////////////////////////////////////////////////////////////////////////
//	QueueNode.java
//
//	Version:           0.3
//
//
//	  authors:
// 	  initials    name                 organization                email
//	---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//
///////////////////////////////////////////////////////////////////////////////////////////////

/**
* Written by Ivelin Georgiev (2004-2006)
* 
*/

/*
	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.
	This library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
	Lesser General Public License for more details.
	
	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
	USA
	
	Contact Info:
		Bruce Donald
		Duke University
		Department of Computer Science
		Levine Science Research Center (LSRC)
		Durham
		NC 27708-0129 
		USA
		brd@cs.duke.edu
	
	If you use or publish any results derived from the use of this
	program please cite:
	Georgiev I, Lilien R, Donald B. "Improved Pruning Algorithms and
	Divide-and-Conquer Strategies for Dead-End Elimination, with Application 
	to Protein Design" Bioinformatics, 22(14): e174-e183, 2006.
	
	Copyright (C) 2006 Ivelin Georgiev, Ryan H. Lilien, and Bruce R. Donald
		
	<signature of Bruce Donald>, 23 Aug, 2006
	Bruce Donald, Professor of Computer Science
*/

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;

public class QueueNode implements Comparable {

	//the f(n) score associated with the current node
	public double fScore;
	
	//corresponding level
	public int level;
	
	//the numbers of the nodes in the considered conformation up to the current level
	public int confSoFar[];
	
	//the number of the corresponding node at that level
	public int nodeNum;
	
	//a pointer to the previous and next nodes in the expansion list
	public QueueNode prevNode;
	public QueueNode nextNode;
	
	//constructor
	QueueNode (int curNode, int curLevel, int curConf[], double fn) {

		nodeNum = curNode;
		level = curLevel;
		
		confSoFar = new int[level+1];
		
		for (int i=0; i<=level; i++){
			confSoFar[i] = curConf[i];
		}	
		fScore = fn;
		
		prevNode = null;
		nextNode = null;	
	}
	//Checks if the two given nodes have the same partially assigned conformation
        private boolean checkConf(QueueNode node1, QueueNode node2){

                if (node1.level!=node2.level) //different level
                        return false;

                for (int l=0; l<node1.confSoFar.length; l++){
                        if (node1.confSoFar[l]!=node2.confSoFar[l])
                                return false;
                }

                //The partially assigned conformations are the same
                return true;
        }

	 public int compareTo(Object otherObject) throws ClassCastException{
               if(!(otherObject instanceof QueueNode)){
                       throw new ClassCastException("A QueueNode object expected.");
               }
               QueueNode other = (QueueNode)otherObject;
               if(this.fScore > other.fScore){
                       return 1;
               }
               else if(this.fScore < other.fScore){
                       return -1;
               }
               else{ // Nodes have the same score
                       if(this.level != other.level || this.nodeNum != other.nodeNum){
                               return 1; // but different levels, this is larger by default
                       }
                       else if(checkConf(this, other)){
                               return 0;
                       }
                       else{ // Two distinct nodes have the same fScore, say this one is larger.
                               return 1;
                       }
               }
       }
}

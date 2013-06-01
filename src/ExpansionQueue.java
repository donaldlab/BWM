///////////////////////////////////////////////////////////////////////////////////////////////
//	ExpansionQueue.java
//
//	Version:           0.3
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
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

/*
 * This queue is ordered in terms of increasing f(n) values of the nodes;
 * 		only the visible nodes are contained in the queue
 */
import java.util.*;

public class ExpansionQueue {
	
	//a pointer to the first node in the expansion list
	public QueueNode curFront;
	
	//number of nodes in the list
	public int numNodes;
	private PriorityQueue<QueueNode> thequeue;
	
	//the unique id of each node in the queue
	public int idNUM;
	
	//constructor
	ExpansionQueue () {
		curFront = null;
		numNodes = 0;
	        thequeue = new PriorityQueue<QueueNode>(50000);
	}
	public void insert(QueueNode newNode){
               thequeue.add(newNode);
               curFront = thequeue.peek();
       }
       public void delete(QueueNode delNode){
               thequeue.remove(delNode);
               curFront = thequeue.peek();
       }
	
	//inserts a new node into the expansion queue
	/*public void insert (QueueNode newNode){
		
		boolean done;
		int i;
		QueueNode curNode = null;
		
		if (numNodes==0){//first node to be inserted
			curFront = newNode;
			numNodes++;
		}
		else {//expansion list is not empty
			curNode = curFront;
			done = false;
			i=0;
			while ((!done)&&(i<numNodes)){
				if (curNode.fScore > newNode.fScore){//insert the new node right before curNode
					
					newNode.nextNode = curNode;
					newNode.prevNode = curNode.prevNode;
					
					if (i!=0)//curNode!=curMinNode
						curNode.prevNode.nextNode = newNode;
					else //we have a new minimum node
						curFront = newNode;
					curNode.prevNode = newNode;					
					
					numNodes++;
					done = true;
				}
				else {
					if (i==numNodes-1){
						newNode.prevNode = curNode;
						curNode.nextNode = newNode;
						numNodes++;
						done = true;
					}
					else {
						curNode = curNode.nextNode;
						i++;
					}
				}
			}
		}
	}

	//Deletes the specified node from the queue: determines which node to delete
	//	based on the fScore, the node level, the node number, and confSoFar[] (only
	//	the last one is sufficient, but it is faster not to check the full length
	//	of confSoFar[] for each node in the queue)
	public void delete (QueueNode delNode){
		
		QueueNode curNode = curFront;
		int i = 0;
		boolean done = false;
		
		while (!done){
			if(delNode.fScore==curNode.fScore){
				if((delNode.level==curNode.level)&&(delNode.nodeNum==curNode.nodeNum)){
					if (checkConf(delNode,curNode)){ //found the node to delete
						if(i!=numNodes-1)
							curNode.nextNode.prevNode = curNode.prevNode;
						if(i!=0)
							curNode.prevNode.nextNode = curNode.nextNode;
						if(i==0)//move the front pointer if we are deleting the min node
							curFront = curNode.nextNode;
						
						numNodes--;
						done = true;
						break;
					}
				}
			}
			if (!done){ //go to the next node in the queue
				curNode = curNode.nextNode;
				i++;
			}
		}
	}*/
	
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
}

package BDAStar;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;

public class BDAStarNode implements Comparable<BDAStarNode> 
{
    
    private PriorityQueue<BDAStarNode> children;
    private Map<Choice, BDAStarNode> childMap;
    private BDAStarNode parent;
    private Choice choice;
    private Position position;
    private Conformation emptyConformation;
    private double score = Double.MAX_VALUE;
    private Conformation prefix;
    
    public BDAStarNode (BDAStarNode parentNode, Choice newChoice, Conformation empty) {
        parent = parentNode;
        choice = newChoice;
        childMap = new HashMap<Choice, BDAStarNode>();
        children = new PriorityQueue<BDAStarNode>();
        emptyConformation = empty;
    }
   

    public void insertConformation(Conformation c, Conformation prefix)
    {
        insertConformation(c, prefix, c.getPositions().toArray(new Position[]{}),0);
    }
    
    private void insertConformation(Conformation c, Conformation prefix, Position[] positions, int index)
    {
        Choice nextChoice = c.getChoiceAt(positions[index]);
        
        if(!childMap.containsKey(nextChoice))
        {
            BDAStarNode newNode = new BDAStarNode(this, nextChoice, emptyConformation);
            newNode.position = positions[index];
            newNode.prefix = prefix;
            newNode.score = c.join(prefix).score();
            childMap.put(nextChoice, newNode);
            children.add(newNode);
            
        	if(index == positions.length -1)
        		heapCheck();
        }
        if(index < positions.length - 1)
            childMap.get(nextChoice).insertConformation(c, prefix, positions, index + 1);
        
    }
    
    public void deleteConformation(Conformation c)
    {
        deleteConformation(c, c.getPositions().toArray(new Position[]{}), 0);
    }
    
    private void deleteConformation(Conformation c, Position[] positions, int index)
    {
        Choice currentChoice = c.getChoiceAt(positions[index]);
        if(!childMap.containsKey(currentChoice))
            return;
        BDAStarNode toDelete = childMap.get(currentChoice);
        deleteConformation(c, positions, index + 1);
        children.remove(toDelete);
    }
    
    public String toString()
    {
    	return super.toString() + " - " + getConformation().toString()+ ":"+getConformation().score();
    }
    
    public Conformation getNextConformation()
    {
        Conformation out = null;
        if(children.size() < 1)
        {
            out = getConformation();
            return out;
        }
        BDAStarNode nextChild = children.poll();
        out = nextChild.getNextConformation();
        if(nextChild.moreConformations())
        	children.add(nextChild);
        else childMap.remove(nextChild.choice);
        return out;
    }
    
    public Conformation peekNextConformation()
    {
        Conformation out = null;
        if(children.size() < 1)
        {
        	return getConformation();
        }
        out = children.peek().peekNextConformation();
        return out;
    }
    
    public Conformation getConformation()
    {
        Conformation out = null;
        if(parent == null)
        {
            out = emptyConformation.copy();
        }
        else 
            out = parent.getConformation();
        if(position != null)
        {
        	out.append(position, choice);
        }
        if(false && children.size() < 1 && prefix != null)
        	out = out.join(prefix);
        return out;
    }
    
    public void heapCheck()
    {
    	if(children.size() < 1) return;
    	for(BDAStarNode child: children)
    	{
    		child.heapCheck();
    	}
    	double lastScore = children.peek().peekNextConformation().score();
    	for(BDAStarNode child: children)
    	{
    		double score = child.nextBestScore();
    		if(score < lastScore)
    		{
    			System.out.println("HEAP OUT OF ORDER: " + score +"," +lastScore);
    			parent.printTree();
    			resort();
    		}
    		lastScore = score;
    	}
    }
    
    public double nextBestScore()
    {
    	//if(children.size() < 1)
    		//return score;
    	return peekNextConformation().score();
    }
    
    public void resort()
    {
    	for(BDAStarNode child: children)
    	{
    		child.resort();
    	}
    	PriorityQueue<BDAStarNode> newChildren = new PriorityQueue<BDAStarNode>();
    	newChildren.addAll(children);
    	children = newChildren;
    }
    
    public void reinsertChain(List<BDAStarNode> path)
    {
        for(int i = 1; i < path.size(); i++)
        {
                BDAStarNode previous = path.get(i);
                previous.children.add(previous.children.poll());
        }
    }
    
    private List<BDAStarNode> peekPartialPath(List<BDAStarNode> path)
    {
        if(path == null)
                path = new LinkedList<BDAStarNode>();
        if(children.size() > 1)
            children.peek().peekPartialPath(path);
        path.add(this);
        return path;
    }
    
    public void printTree()
    {
    	printTree("");
    }
    
    public void printTree(String prefix)
    {
    	Conformation peeked = peekNextConformation();
    	String out = prefix + peeked +" : "+peeked.score();
    	System.out.println(out);
    	for(BDAStarNode child : children)
    	{
    		child.printTree(prefix + "+--");
    	}
    }

    @Override
    public int compareTo (BDAStarNode arg0) {
        double diff = nextBestScore() - arg0.nextBestScore();
        if(diff < 0) return -1;
        if(diff > 0) return 1;
        return 0;
    }

    public boolean moreConformations () {
        return children.size() > 0;
    }


}

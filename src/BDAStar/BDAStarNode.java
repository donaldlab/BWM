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
    private double score;
    
    public BDAStarNode (BDAStarNode parentNode, Choice newChoice, Conformation empty) {
        parent = parentNode;
        choice = newChoice;
        childMap = new HashMap<Choice, BDAStarNode>();
        children = new PriorityQueue<BDAStarNode>();
        emptyConformation = empty;
    }

    public void insertConformation(Conformation c)
    {
        insertConformation(c, c.getPositions().toArray(new Position[]{}),0);
    }
    
    private void insertConformation(Conformation c, Position[] positions, int index)
    {
        Choice nextChoice = c.getChoiceAt(positions[index]);
        
        if(!childMap.containsKey(nextChoice))
        {
            BDAStarNode newNode = new BDAStarNode(this, nextChoice, emptyConformation);
            children.add(newNode);
            newNode.position = positions[index];
            childMap.put(nextChoice, newNode);
            if(index + 1 == positions.length)
                newNode.score = c.score();
        }
        if(index < positions.length - 1)
            childMap.get(nextChoice).insertConformation(c, positions, index + 1);
        score = children.peek().getConformation().score();
        
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
    
    public Conformation getNextConformation()
    {
        Conformation out = null;
        if(children.size() < 1)
        {
            out = emptyConformation;
            if(position!= null)
            out.append(position, choice);
            out.assignScore(score);
            return out;
        }
        BDAStarNode nextChild = children.poll();
        out = nextChild.getNextConformation();
        if(position != null)
        out.append(position, choice);
        if(nextChild.moreConformations())
        	children.add(nextChild);
        return out;
    }
    
    public Conformation peekNextConformation()
    {
        Conformation out = null;
        if(children.size() < 1)
        {
            out = emptyConformation;
            if(position != null)
            out.append(position, choice);
            return out;
        }
        out = children.peek().peekNextConformation();
        if(position != null)
        out.append(position, choice);
        return out;
    }
    
    public Conformation getConformation()
    {
        Conformation out = null;
        if(parent == null)
        {
            out = emptyConformation;
        }
        else 
            out = parent.getConformation();
        if(position != null)
        out.append(position, choice);
        return out;
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

    @Override
    public int compareTo (BDAStarNode arg0) {
        double diff = getConformation().score() - arg0.getConformation().score();
        if(diff < 0) return -1;
        if(diff > 0) return 1;
        return 0;
    }

    public boolean moreConformations () {
        return children.size() > 0;
    }


}

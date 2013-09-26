package BDAStar;

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
    private BWMSolutionSpace space;
    private double score;
    
    public BDAStarNode (BDAStarNode parentNode, Choice newChoice) {
        parent = parentNode;
        choice = newChoice;
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
            BDAStarNode newNode = new BDAStarNode(this, nextChoice);
            children.add(newNode);
            childMap.put(nextChoice, newNode);
            if(index + 1 == positions.length)
                newNode.score = c.score();
        }
        score = children.peek().getConformation().score();
        childMap.get(c.getChoiceAt(positions[index])).insertConformation(c, positions, index + 1);
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
            out = space.getEmptyConformation();
            ProteinConformation conf = (ProteinConformation) out;
            out.append(position, choice);
            conf.assignScore(score);
            return out;
        }
        out = children.poll().getNextConformation();
        out.append(position, choice);
        return out;
    }
    
    public Conformation peekNextConformation()
    {
        Conformation out = null;
        if(children.size() < 1)
        {
            out = space.getEmptyConformation();
            out.append(position, choice);
            return out;
        }
        out = children.peek().getNextConformation();
        out.append(position, choice);
        return out;
    }
    
    public Conformation getConformation()
    {
        Conformation out = null;
        if(parent == null)
        {
            out = space.getEmptyConformation();
        }
        else 
            out = parent.getConformation();
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

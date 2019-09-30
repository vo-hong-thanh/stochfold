/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package model.secondarystructure;

import constants.RINGNODE_CODE;

/**
 *
 * @author vot2
 */
public class RingNode {
    private int positionInSequence;
    
    private char type;
    
    private RingNode next;
    private RingNode prev;
    private RingNode up;
    private RingNode down;
    
    public RingNode(){
        positionInSequence = -1;
        type = RINGNODE_CODE.UNKNOWN;
        next = null;
        prev = null;
        up = null;
        down = null;
    }
    
    //position
    public void setPosition(int _position)
    {
        positionInSequence = _position;
    }
    
    public int getPosition()
    {
        return positionInSequence;
    }
    
    //type of node
    public void setType(char _type)
    {
        type = _type;
    }
    
    public char getType()
    {
        return type;
    }
    
    //next
    public void setNext(RingNode _next)
    {
        next = _next;
    }
    
    public RingNode getNext()
    {
        return next;
    }
    
    //prev
    public void setPrev(RingNode _prev)
    {
        prev = _prev;
    }
    
    public RingNode getPrev()
    {
        return prev;
    }
    
    //up
    public void setUp(RingNode _up)
    {
        up = _up;
    }
    
    public RingNode getUp()
    {
        return up;
    }
    
    //down
    public void setDown(RingNode _down)
    {
        down = _down;
    }
    
    public RingNode getDown()
    {
        return down;
    }
    
    @Override
    public String toString(){
        StringBuilder s = new StringBuilder("RingNode (position = " + positionInSequence +")");
        s.append("\n");
        s.append(" - Pairing type = " + type);
        s.append("\n");
        s.append(" - Links to (");
        s.append("next " + (next == null? "null" : "at position: " + next.getPosition()) );
        s.append(", previous " + (prev == null ? "null" : "at position: " + prev.getPosition()) );
        s.append(", down " + (down == null? "null" : "at position: " + down.getPosition()) );
        s.append(", up " + (up == null? "null" : "at position: " + up.getPosition()) );
        s.append(")");
        return s.toString();
    }
    
}

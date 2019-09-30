/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author vot2
 */
public class Cache<K, V> extends LinkedHashMap<K, V>{
    private static final float loadFactor = 0.75f;
    private static final boolean accessOrder = true;
    
    private final int cacheSize;
    
    public Cache(int _cacheSize) {
        super(_cacheSize, loadFactor, accessOrder);
        cacheSize = _cacheSize;
    }

    @Override
    protected boolean removeEldestEntry(Map.Entry<K, V> eldest) {
        return size() > cacheSize;
    }
}

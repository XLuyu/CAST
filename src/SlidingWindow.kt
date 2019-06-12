package edu.nus

class SlidingWindowArray(val size:Int){
    private val window = Array(size){ arrayOf(0,0,0,0,0) }
    fun inc(i:Int,j:Int) { window[i%size][j] += 1 }
    fun getAndClean(i:Int): Array<Int> {
        val s = window[i%size]
        window[i%size] = arrayOf(0,0,0,0,0)
        return s
    }
}

class SlidingWindowInt(val size:Int){
    private val window = Array(size){0}
    fun inc(i:Int) { window[i%size] += 1 }
    fun getAndClean(i:Int): Int {
        val s = window[i%size]
        window[i%size] = 0
        return s
    }
}
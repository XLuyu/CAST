package edu.nus

import kotlinx.coroutines.*
import me.tongfei.progressbar.ProgressBar

class BamFileDetector(filename:String): BamFileScanner(filename){
    private var coverageStat = mutableMapOf<Int,Long>()

    override fun get(cid:Int,pos:Int): Genotype {
        val genotype = super.get(cid, pos)
        if (genotype.isReliable) coverageStat[genotype.sum] =  coverageStat.getOrDefault(genotype.sum,0) + 1
        return genotype
    }
    fun getCoverage(): Double {
        val lowerbound = coverageStat.map{it.key*it.value}.sum()/coverageStat.values.sum().toDouble() //TODO: too high, use Int
        return coverageStat.filter{it.key>=lowerbound}.maxBy{it.value}?.key?.toDouble() ?:0.0
    }
}
class DepthDetector(filenames:List<String>){
    private val bamScanners = filenames.map{ BamFileDetector(it) }
    private val headers = bamScanners.map{ it.header }
    private val header = headers[0]
    private var cid = header.sequences.maxBy { it.sequenceLength }!!.sequenceIndex
    private val maxpos = header.getSequence(cid).sequenceLength
    private val pb = ProgressBar("Pre-check", maxpos.toLong())
    init { if (headers.any {it != header}) throw Exception("[Error] Headers in input files are inconsistent") }

    val getCoverage = {
        for (pos in 1..maxpos){
            if (pos==1 || pos>bamScanners[0].cached)
                runBlocking {
                    bamScanners.map { async (Dispatchers.Default) { it.get(cid, pos) } }.map { it.await() }
                }
            else
                bamScanners.map { it.get(cid, pos)}
            if (pos%1000==0) pb.stepBy(1000)
        }
        pb.close()
        bamScanners.map{ it.getCoverage() }
    }()
}

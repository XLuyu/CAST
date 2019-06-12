package edu.nus

import htsjdk.samtools.*
import htsjdk.samtools.SAMRecord
import kotlinx.coroutines.*
import me.tongfei.progressbar.ProgressBar
import kotlin.system.exitProcess

class BamFileDetector(val filename:String): BamFileScanner(filename){
    var coverageStat = mutableMapOf<Int,Long>()

    override fun get(cid:Int,pos:Int): Genotype {
        val genotype = super<BamFileScanner>.get(cid, pos)
        if (genotype.isReliable) coverageStat[genotype.sum] =  coverageStat.getOrDefault(genotype.sum,0) + 1
        return genotype
    }
    fun getCoverage(): Double {
        val lowerbound = coverageStat.map{it.key*it.value}.sum()/coverageStat.values.sum().toDouble() //TODO: too high, use Int
        return coverageStat.filter{it.key>=lowerbound}.maxBy{it.value}?.key?.toDouble() ?:0.0
    }
}
class DepthDetector(filenames:List<String>){
    val bamScanners = filenames.map{ BamFileDetector(it) }
    val headers = bamScanners.map{ it.header }
    val header = headers[0]
    val contigNum = header.sequences.size
    var cid = (0 until contigNum).maxBy{header.getSequence(it).sequenceLength }!!
    var pos = 0
    private val maxpos = header.getSequence(cid).sequenceLength
    val pb = ProgressBar("Pre-check", maxpos.toLong())
    init { if (headers.any {it != header}) throw Exception("[Error] Headers in input files are inconsistent") }
    val getCoverage = {
        while (pos<maxpos){
            pos += 1
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

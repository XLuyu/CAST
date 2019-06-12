package edu.nus

import htsjdk.samtools.*
import me.tongfei.progressbar.ProgressBar
import java.io.File
import kotlinx.coroutines.*

open class BamFileScanner(filename:String){
    val bamFile = SamReaderFactory.makeDefault().open(File(filename))
    val header = bamFile.fileHeader.sequenceDictionary
    val fileIterator = bamFile.iterator()
    var currentRead:SAMRecord? = fileIterator.next()
    var sw = SlidingWindowArray(100000)
    var badRegion = SlidingWindowInt(100000)

    val code = mapOf<Char,Int>('A' to 0,'C' to 1,'G' to 2,'T' to 3)
    var cached = 0
    var cigarRegex = Regex("""(\d+)[SHDI]""")

    private fun isBadRead(read:SAMRecord) = read.mappingQuality <30 || read.mateUnmappedFlag || read.mateReferenceName !=read.referenceName ||
        cigarRegex.findAll(read.cigarString).map{it.groups[1]!!.value.toInt() }.map{it*it-1}.sum()+read.getIntegerAttribute("NM")>0.1*read.readLength ||
        read.hasAttribute("XA")
    fun updateSpanReadsByPosition(cid:Int,pos:Int): Unit {
        cached = pos + sw.size/2
        while (currentRead!=null) {
            if (currentRead!!.readUnmappedFlag) { currentRead = fileIterator.next(); continue }
            val icid = currentRead!!.referenceIndex
            val ipos = currentRead!!.alignmentStart
            if (icid < cid) { currentRead = fileIterator.next(); continue }
            if (icid == cid && ipos <= cached) {
                val record = currentRead!!
                currentRead = fileIterator.next();
                val bad = isBadRead(record)
                val read = record.readString
                val qual = record.baseQualities
                for (i in record.start..record.end) {
                    val posInRead = record.getReadPositionAtReferencePosition(i)
                    if (posInRead == 0) {
                        sw.inc(i, 4)
                    } else {
                        val ch = read[posInRead - 1]
                        if (code.contains(ch) && qual[posInRead - 1] >= 20) sw.inc(i, code.getValue(ch))
                    }
                    if (bad) badRegion.inc(i)
                }
            } else break
        }
    }
    open fun get(cid:Int, pos:Int): Genotype {
        if (pos==1) cached = 0
        if (pos>cached) updateSpanReadsByPosition(cid,pos)
        val badCount = badRegion.getAndClean(pos)
        val count = sw.getAndClean(pos)
        return Genotype(count.map{ it.toDouble() },badCount)
    }
}
class BamFilesParallelScanner(filenames:List<String>){
    val bamScanners = filenames.map{ BamFileScanner(it) }
    val headers = bamScanners.map{ it.header }
    val header = headers[0]
    var cid = 0
    var pos = 0
    private var maxpos = header.getSequence(cid).sequenceLength
    val pb = ProgressBar("Scan Bam files", header.referenceLength)
    init { nextPosition() }
    fun nextPosition(): Boolean { // return true if on same contig
        pos += 1
        Dispatchers.Default
        if (pos>maxpos){
            pb.stepBy(pos%1000L).extraMessage = "finished ${getChrom()}"
            pos = 1
            do {
                cid += 1
                if (cid<header.size()) {
                    maxpos = header.getSequence(cid).sequenceLength
                    if (maxpos<5000) pb.stepBy(maxpos.toLong())
                }
                else pb.close()
            } while (cid<header.size() && maxpos<5000)
        }
        if (pos%1000==0) pb.stepBy(1000)
        return pos!=1
    }
    fun get() =
        if (pos==1 || pos>bamScanners[0].cached)
            runBlocking {
                bamScanners.map { async (Dispatchers.Default) { it.get(cid, pos) } }.map { it.await() }
            }
        else
            bamScanners.map { it.get(cid, pos)}
    fun getChrom() = header.getSequence(cid).sequenceName
    fun hasNext() = cid<header.size()
}

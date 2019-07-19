package edu.nus

import com.github.ajalt.clikt.core.*
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.arguments.multiple
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.file
import org.tc33.jheatchart.HeatChart
import java.io.File

class CAST(private val args: Array<String>) : CliktCommand() {
    private val outDir by option("-o", help="output directory [CAST_output]").file().default(File("CAST_output"))
    private val heatmap by option(help="output distance matrix heatmap [None]").flag()
    private val draftFile by argument(help="draft assembly (usually parent) to be improved").file(true)
    private val bamFiles by argument(help=".bam files (usually progeny) to provide genetic information, one file for each individual")
        .multiple(true)

    override fun run() {
        outDir.mkdirs()
        File(outDir,"cmd.log").writeText(args.joinToString(" ","cmd: ","\n"))
        val genotyper = DistanceGenotyper(bamFiles, outDir, heatmap)
        val reportFile = genotyper.run()
        val stringer = Scaffolder(draftFile, reportFile, outDir)
        stringer.correctAndScaffoldFasta()
    }
}
fun main(args: Array<String>) = CAST(args).main(if (args.isEmpty()) arrayOf("--help") else args)
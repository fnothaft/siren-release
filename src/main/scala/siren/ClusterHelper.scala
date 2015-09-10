package siren

import java.io._

object ClusterHelper {
  def main(args: Array[String]) : Unit = {
    args match {
      case Array("prepareClusters", clusterFilename, outputFilename, genomeFilename) => {
        // Open cluster file
        val source = scala.io.Source.fromFile(clusterFilename)
        val lines = source.getLines
        
        // Load genome
        val genome = FASTA.read(genomeFilename)

        // Open output file
        val fw = new java.io.FileWriter(outputFilename)
        val bw = new java.io.BufferedWriter(fw)

        var i = 0
        lines.foreach(l => {
          if (i % 1000 == 0)
            println("Parsing line " + i + "...")

          // Skip first line
          if (i > 0) {
            val entries = l.split(" ")
            val numEntries = entries.length
            // entries(0) is cluster size
            val clusterId = entries(1).toLong
            val (idPiece, idPos) = genome.getLocation(clusterId)
            (1 to numEntries).foreach(n => {
              val (piece, pos) = genome.getLocation(entries(n).toLong)
              bw.write(List(piece, pos, idPiece, idPos).mkString("\t"))
              bw.newLine
            })

            i += 1
          }
        })
        
        // Close output file
        bw.close
      }
      case Array("countErrors", clusterFilename, samFilename, genomeFilename) => {
        // load cluster file & create hash map of cluster positions to cluster id
        /*
        println("Loading cluster file...")
        val source = scala.io.Source.fromFile(clusterFilename)
        val lines = source.getLines
        val clusterMembership = scala.collection.mutable.Map[(String, Long), (String, Long)]()

        var i = 0
        lines.foreach(l => {
          if (i % 10000000 == 0)
            println("Parsing line " + i + "...")

          val entries = l.split("\t")
          assert(entries.length == 4)
          val pos = (entries(0), entries(1).toLong)
          val id = (entries(2), entries(3).toLong)
          val returnVal = clusterMembership.getOrElseUpdate(pos, id)
          if (returnVal != id) {
            println("Warning:  pos " + pos + " was already in clusterMembership map.")
            println("Old cluster: " + returnVal + ", new cluster: " + id)
          }

          i += 1
        })
        println("Done loading cluster file.")
        */
        // Load genome
        println("Loading genome...")
        val genome = FASTA.read(genomeFilename)

        // Load clusters (will convert to 1-indexed)
        val clusterMembership = getClusterMembers(clusterFilename)
        
        // iterate over sam file
        println("Traversing sam file...")
        
        // for each entry
        // is it a miss?
        // yes => is it part of a cluster?
        val source2 = scala.io.Source.fromFile(samFilename)
        val lines2 = source2.getLines
        
        var numErrors = 0
        var numErrorsInClusters = 0
        
        val numErrorsByMapq = Array.fill(71)(0) // to accommodate 0 to 70
        val numErrorsInClustersByMapq = Array.fill(71)(0)
        
        var i = 0
        lines2.foreach(l => {
          if (i % 10000 == 0)
            println("Parsing line " + i + "...")
          
          val entries = l.split("\t")
          entries(12) match {
            case "mis" => {
              val mapq = entries(4).toInt
              numErrors += 1
              numErrorsByMapq(mapq) += 1
              
              val piece = entries(2)
              val offset = entries(3).toLong
              /*
              clusterMembership.get((piece, offset)) match {
                case Some((idPiece, idOffset)) => {
                  numErrorsInClusters += 1
                  numErrorsInClustersByMapq(mapq) += 1
                }
                case None => 
              }
              */
              val absPos = genome.getAbsPos(piece, offset)
              if (clusterMembership.contains(absPos)) {
                numErrorsInClusters += 1
                numErrorsInClustersByMapq(mapq) += 1
              }
            }
            case _ =>
          }
            
          i += 1
        })
        
        println("Done traversing sam file.")
        
        // print out # entries
        // print out # errors
        // print out # errors that were in clusters
        println("# entries in SAM file: " + i)
        printf("# errors: %d (%.3f%%)\n", numErrors, (numErrors.toDouble / i * 100.0))
        if (numErrors > 0)
          printf("# errors in clusters: %d (%.3f%%)\n", numErrorsInClusters, (numErrorsInClusters.toDouble / numErrors * 100.0))
          
        // print # errors, # errors in clusters (cumulative) vs. mapq
        var numErrorsSoFar = 0
        var numErrorsInClustersSoFar = 0
        
        println("# errors, # errors in clusters (cumulative) vs. mapq")
        println(List("MAPQ", "N Errors", "N Errors in Clusters").mkString("\t"))
        
        var mapq = 70
        while (mapq >= 0) {
          numErrorsSoFar += numErrorsByMapq(mapq)
          numErrorsInClustersSoFar += numErrorsInClustersByMapq(mapq)
          println(List(mapq, numErrorsSoFar, numErrorsInClustersSoFar).mkString("\t"))
          
          mapq -= 1
        }
      }
      case Array("countRepeatMaskerPositions", repeatMaskerFilename, readLenStr, genomeFilename) => {
        val readLen = readLenStr.toInt
        
        // load repeat masker filename
        val source = scala.io.Source.fromFile(repeatMaskerFilename)
        val lines = source.getLines

        var numSubstrings = 0L

        var i = 0
        lines.foreach(l => {
          if (i % 1000000 == 0)
            println("Parsing line " + i + "...")
            
          // Skip header (first 3 lines)
          if (i > 2) {
            // parse line
            
            /*
            scala> " 1504   1.3  0.4  1.3  chr1        10001   10468 (249240153) +  (CCCTAA)n      Simple_repeat            1  463    (0)      1"
            res0: java.lang.String = " 1504   1.3  0.4  1.3  chr1        10001   10468 (249240153) +  (CCCTAA)n      Simple_repeat            1  463    (0)      1"

            scala> res0.trim
            res4: java.lang.String = 1504   1.3  0.4  1.3  chr1        10001   10468 (249240153) +  (CCCTAA)n      Simple_repeat            1  463    (0)      1

            scala> res4.split(" +")
            res7: Array[java.lang.String] = Array(1504, 1.3, 0.4, 1.3, chr1, 10001, 10468, (249240153), +, (CCCTAA)n, Simple_repeat, 1, 463, (0), 1)
            */

            val entries = l.trim.split(" +")
            val highEnd = entries(6).toLong
            val lowEnd = entries(5).toLong
            
            val numSubstringsInRange = (highEnd - lowEnd + 1) - readLen + 1
            numSubstrings += numSubstringsInRange
          }
          
          i += 1
        })

        // Load genome
        val genome = FASTA.read(genomeFilename)
        val genomeSize = genome.totalSize

        // print stats
        println("# substrings in repeats: " + numSubstrings)
        println("# positions in genome: " + genomeSize)
        printf("%% of genome in repeats: %.3f%%\n", (numSubstrings.toDouble / genomeSize * 100.0))
      }
      case Array("makeRocCurve", samFilename) => {
        val source = scala.io.Source.fromFile(samFilename)
        val lines = source.getLines
        
        var numCorrect = 0
        var numIncorrect = 0
        var numUnaligned = 0
        
        val numCorrectByMapq = Array.fill(71)(0) // to accommodate 0 to 70
        val numIncorrectByMapq = Array.fill(71)(0)
        
        var i = 0
        lines.foreach(l => {
          if (i % 10000 == 0)
            println("Parsing line " + i + "...")
          
          val entries = l.split("\t")
          val mapq = entries(4).toInt
          entries(12) match {
            case "mis" => {
              numIncorrect += 1
              numIncorrectByMapq(mapq) += 1
            }
            case "un" => {
              numUnaligned += 1
            }
            case _ => {
              numCorrect += 1
              numCorrectByMapq(mapq) += 1
            }
          }
          
          i += 1
        })
        
        // Print out ROC curve
        var numCorrectSoFar = 0
        var numIncorrectSoFar = 0
        
        println(List("MAPQ", "numCorrect", "numIncorrect").mkString("\t"))
        
        var mapq = 70
        while (mapq >= 0) {
          numCorrectSoFar += numCorrectByMapq(mapq)
          numIncorrectSoFar += numIncorrectByMapq(mapq)
          println(List(mapq, numCorrectSoFar, numIncorrectSoFar).mkString("\t"))
          
          mapq -= 1
        }
        
      }
      case Array("outputReadsInClusters", clusterFilename, samFilename, outputSamFilename, genomeFilename) => {
        // Read in clusters; add to set
        println("Loading cluster file...")
        val source = scala.io.Source.fromFile(clusterFilename)
        val lines = source.getLines
        val clusterMembership = scala.collection.mutable.Set[Long]()

        var i = 0
        lines.foreach(l => {
          if (i % 10000000 == 0)
            println("Parsing line " + i + "...")

          // Skip first line
          if (i > 0) {
            val entries = l.split(" ")
            // entries(0) is cluster size
            (1 until entries.length).foreach(n => clusterMembership += (entries(n).toLong + 1 /* because 0-indexed; want 1-indexed */ ))
          }
            
          i += 1
        })
        println("Done loading cluster file.")
        
        // Open output file
        val fw = new java.io.FileWriter(outputSamFilename)
        val bw = new java.io.BufferedWriter(fw)
        
        // Load genome
        val genome = FASTA.read(genomeFilename)

        // Iterate over sam file
        // Is a read's aligned pos in a cluster?
        // Yes => output to file
        val source2 = scala.io.Source.fromFile(samFilename)
        val lines2 = source2.getLines
        
        var numReadsInClusters = 0
        
        i = 0
        lines2.foreach(l => {
          if (i % 20000 == 0)
            println("Parsing line " + i + "...")
          
          val entries = l.split("\t")
          val piece = entries(2)
          val offset = entries(3).toLong
          val absPos = genome.getAbsPos(piece, offset)
          if (clusterMembership.contains(absPos)) {
            numReadsInClusters += 1
            bw.write(l)
            bw.newLine
          }
            
          i += 1
        })
        
        println("Done traversing sam file.")
        
        // Close output file
        bw.close
        
        // Print # reads in clusters
        printf("Number of reads in clusters: %d (%.3f%%)\n", numReadsInClusters, (numReadsInClusters.toDouble / i * 100.0))
      }
      case Array("compareSnapAndNovo", clusterFilename, snapSam, novoSam, genomeFilename) => {
        // Load genome
        println("Loading genome...")
        val genome = FASTA.read(genomeFilename)
        
        // Load clusters
        println("Loading clusters...")
        val clusterMembership = getClusterMembers(clusterFilename)
        
        // Iterate over snap, novo sam files
        val snapCorrect = Array.fill(71)(0) // to accommodate 0 to 70
        val snapIncorrect = Array.fill(71)(0)
        val novoCorrect = Array.fill(71)(0)
        val novoIncorrect = Array.fill(71)(0)
        val hybridCorrect = Array.fill(71)(0)
        val hybridIncorrect = Array.fill(71)(0)
        
        val readIdsInClusters = scala.collection.mutable.Set[String]()
        var numReadsForHybrid = 0
        
        // populate snap roc & non-cluster part of hybrid roc
        println("Loading snap sam file...")

        val source = scala.io.Source.fromFile(snapSam)
        val lines = source.getLines

        var i = 0
        lines.foreach(l => {
          if (i % 20000 == 0)
            println("Parsing line " + i + "...")
          val entries = l.split("\t")

          // record entry in snap roc
          val mapq = entries(4).toInt
          entries(12) match {
            case "mis" => snapIncorrect(mapq) += 1
            case _ => snapCorrect(mapq) += 1
          }

          // check to see if this entry is in a cluster
          // yes => skip
          // no => record in hybrid roc (ie, trust snap whenever read is NOT in cluster)
          val piece = entries(2)
          val offset = entries(3).toLong
          val absPos = genome.getAbsPos(piece, offset)
	        val id = entries(0)          

          if (clusterMembership.contains(absPos)) {
            readIdsInClusters += id
          } else {
            numReadsForHybrid += 1
            entries(12) match {
              case "mis" => hybridIncorrect(mapq) += 1
              case _ => hybridCorrect(mapq) += 1
            }
          }

          i += 1
        })
        
        println("Done traversing snap sam file.")
        println("Loading novo sam file...")

        val source2 = scala.io.Source.fromFile(novoSam)
        val lines2 = source2.getLines

        i = 0
        lines2.foreach(l => {
          if (i % 20000 == 0)
            println("Parsing line " + i + "...")
          val entries = l.split("\t")

          // record entry in novo roc
          val mapq = entries(4).toInt
          entries(12) match {
            case "mis" => novoIncorrect(mapq) += 1
            case _ => novoCorrect(mapq) += 1
          }

          // check to see if snap skipped this entry
          // yes => record in hybrid roc (ie, trust novo whenever read IS in cluster)
          // no => skip
          val id = entries(0)
          if (readIdsInClusters.contains(id)) {
            numReadsForHybrid += 1
            entries(12) match {
              case "mis" => hybridIncorrect(mapq) += 1
              case _ => hybridCorrect(mapq) += 1
            }
          }

          i += 1
        })
        println("Done traversing novo sam file.")
        
        // Print out ROC curves for snap, novo, & hybrid
        println("# reads for hybrid: " + numReadsForHybrid)
        
        var snapCorrectSoFar = 0
        var snapIncorrectSoFar = 0
        var novoCorrectSoFar = 0
        var novoIncorrectSoFar = 0
        var hybridCorrectSoFar = 0
        var hybridIncorrectSoFar = 0
        
        println(List("", "SNAP", "", "Novoalign", "", "Hybrid", "").mkString("\t"))
        println(List("MAPQ", "numCorrect", "numIncorrect", "numCorrect", "numIncorrect", "numCorrect", "numIncorrect").mkString("\t"))
        
        var mapq = 70
        while (mapq >= 0) {
          snapCorrectSoFar += snapCorrect(mapq)
          snapIncorrectSoFar += snapIncorrect(mapq)
          
          novoCorrectSoFar += novoCorrect(mapq)
          novoIncorrectSoFar += novoIncorrect(mapq)

          hybridCorrectSoFar += hybridCorrect(mapq)
          hybridIncorrectSoFar += hybridIncorrect(mapq)
          
          println(List(mapq, snapCorrectSoFar, snapIncorrectSoFar, novoCorrectSoFar, novoIncorrectSoFar, hybridCorrectSoFar, hybridIncorrectSoFar).mkString("\t"))
          
          mapq -= 1
        }
        
        
      }
      case Array("summarizeReadErrors", samFile, outFile) => {
        // Open output file
        val fw = new java.io.FileWriter(outFile)
        val bw = new java.io.BufferedWriter(fw)

        // iterate over sam file
        // print out (id \t -1/0/1); -1 if error, 0 if unaligned, 1 if correct
        println("Loading sam file...")

        val source = scala.io.Source.fromFile(samFile)
        val lines = source.getLines

        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)

          // check to see if it was aligned correctly
          val code = 
          entries(12) match {
            case "mis" => -1
            case "un" => 0
            case _ => 1
          }
          
          bw.write(id + "\t" + code.toString)
          bw.newLine
        })
        
        bw.close
      }
      case Array("summarizeReadClusterMembership", clusterFile, samFile, genomeFile, outFile) => {
        // Open output file
        val fw = new java.io.FileWriter(outFile)
        val bw = new java.io.BufferedWriter(fw)

        // Load genome
        val genome = FASTA.read(genomeFile)
        
        // Load clusters
        println("Loading clusters...")
        val clusterMembership = getClusterMembers(clusterFile)
        
        // Iterate over sam file
        // print out (id \t 0/1); 0 if not in cluster, 1 if in cluster
        val source = scala.io.Source.fromFile(samFile)
        val lines = source.getLines

        lines.foreach(l => {
	        if (l(0) != '@') {
            val entries = l.split("\t")
            val id = entries(0)
          
            // check to see if this entry is in a cluster
            val piece = entries(2)
            val offset = entries(3).toLong
            val absPos = genome.getAbsPos(piece, offset)
          
            val code = 
            if (clusterMembership.contains(absPos)) 1
            else 0
          
            bw.write(id + "\t" + code.toString)
            bw.newLine
	        }
        })
        
        bw.close
      }
      case Array("analyzeAligners", alignerFile, codeStr) => {
        val code = codeStr.toInt

        // load in aligner file
        val source = scala.io.Source.fromFile(alignerFile)
        val lines = source.getLines

        // look for reads that match the code of interest
        // eg, if code is 1, look for reads where 3 of the aligners report code 1 for that read
        val alignersToCodeInClusters = Array.fill(5)(0) // 0 => no aligners ... 4 => all aligners
        val alignersToCodeNotInClusters = Array.fill(5)(0)
        
        var i = 0
        lines.foreach(l => {
          // skip first line
          if (i > 0) {
            val entries = l.split("\t")
            
            // first column is read id; next 4 columns are codes for the 4 aligners; last column is for cluster membership
            assert(entries.length == 6)
            val numToCode = entries.toList.slice(1, entries.length - 1).map(_.toInt).count(_ == code)

            if (entries(entries.length - 1).toInt == 0) // not in cluster
              alignersToCodeNotInClusters(numToCode) += 1
            else  // in cluster
              alignersToCodeInClusters(numToCode) += 1
          }

          i += 1
        })

        // print summary
        val inClustersSum = alignersToCodeInClusters.sum
        val notInClustersSum = alignersToCodeNotInClusters.sum
        println("Total (in clusters): " + inClustersSum)
        println("Total (not in clusters): " + notInClustersSum)
        println("Total (overall): " + (inClustersSum + notInClustersSum))
        println("NumAligners\tToCodeInClusters\tToCodeNotInClusters\tPercentToCodeInClusters")
        (0 until alignersToCodeInClusters.length).foreach(i => {
          //println("Reads where " + i + " aligners are to code (in clusters): " + alignersToCodeInClusters(i))
          //println("Reads where " + i + " aligners are to code (not in clusters): " + alignersToCodeNotInClusters(i))
          val inClusters = alignersToCodeInClusters(i)
          val notInClusters = alignersToCodeNotInClusters(i)
          println(i + "\t" + inClusters + "\t" + notInClusters + "\t" + (inClusters.toDouble / (inClusters + notInClusters) * 100.0))
        })
      }
      case Array("splitSam", numReadsStr, snapSam, novoSam, outPrefix) => {
        val numReads = numReadsStr.toInt
        
        // for each read, want to know snap & novo result (either correct: 1, un: 0, mis: -1)
        val snapCodes = Array.fill(numReads)(0)
        val novoCodes = Array.fill(numReads)(0)

        // read in snap sam
        println("Reading in snap sam...")
        var source = scala.io.Source.fromFile(snapSam)
        var lines = source.getLines

        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val idAsInt = masonIdToInt(id)
          
          val code = entries(12) match {
            case "mis" => -1
            case "un" => 0
            case _ => 1
          }
          
          snapCodes(idAsInt) = code
        })

        println("Reading in novo sam...")
        source = scala.io.Source.fromFile(novoSam)
        lines = source.getLines
        
        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val idAsInt = masonIdToInt(id)
          
          val code = entries(12) match {
            case "mis" => -1
            case "un" => 0
            case _ => 1
          }
          
          novoCodes(idAsInt) = code
        })
        
        // then, for each read, decide whether to direct towards snap file or novo file
        val snapOut = FASTQ.writer(outPrefix + ".snap.fastq")
        val novoOut = FASTQ.writer(outPrefix + ".novo.fastq")
        
        // reopen snap sam
        source = scala.io.Source.fromFile(snapSam)
        lines = source.getLines
        
        var readsToSnap = 0
        var readsToNovo = 0
        
        val caseCounts = Array.fill(9)(0)
        
        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val seq = entries(9)
          val quality = entries(10)
          val read = new Read(id.getBytes, seq.getBytes, quality.getBytes)

          val i = masonIdToInt(id)
          val snapCode = snapCodes(i)
          val novoCode = novoCodes(i)
          
          // 0 => snap; 1 => novo
          /*
          val split = (snapCode, novoCode) match {
            case (1, _) => 0
            case (0, 1) => 1
            case (0, 0) => 0
            case (0, -1) => 0
            case (-1, 1) => 1
            case (-1, 0) => 1
            case (-1, -1) => 0
          }
          */
          val split = (snapCode, novoCode) match {
            case (1, _) => {
              novoCode match {
                case 1 => caseCounts(0) += 1
                case 0 => caseCounts(1) += 1
                case -1 => caseCounts(2) += 1
              }
              0
            }
            case (0, 1) => {
              caseCounts(3) += 1
              1
            }
            case (0, 0) => {
              caseCounts(4) += 1
              0
            }
            case (0, -1) => {
              caseCounts(5) += 1
              0
            }
            case (-1, 1) => {
              caseCounts(6) += 1
              1
            }
            case (-1, 0) => {
              caseCounts(7) += 1
              1
            }
            case (-1, -1) => {
              caseCounts(8) += 1
              0
            }
          }

          // write out to appropriate fastq file
          if (split == 0) {
            snapOut.write(read)
            readsToSnap += 1
          }
          else {
            novoOut.write(read)
            readsToNovo += 1
          }
        })
        
        snapOut.close
        novoOut.close
        
        // print summary
        println("Printed " + readsToSnap + " reads to snap fastq.")
        println("Printed " + readsToNovo + " reads to novo fastq.")
        println("Printed " + (readsToSnap + readsToNovo) + " total reads to fastq.")
        println("Original # reads: " + numReads)
        
        println("Case counts: " + caseCounts.mkString(", "))
        println("Sum (case counts): " + caseCounts.sum)
      }
      case Array("outputErrorsNotInClusters", clusterFilename, samFilename, outputSamFilename, genomeFilename) => {
        // open output file
        val fw = new java.io.FileWriter(outputSamFilename)
        val bw = new java.io.BufferedWriter(fw)

        // load genome
        println("Loading clusters...")
        val genome = FASTA.read(genomeFilename)

        // load clusters
        println("Loading clusters...")
        val clusterMembership = getClusterMembers(clusterFilename)

        // read in sam file
        val source = scala.io.Source.fromFile(samFilename)
        val lines = source.getLines
        
        var numErrors = 0
        var numErrorsNotInClusters = 0
        
        // for each read
        println("Loading sam file...")
        lines.foreach(l => {
          val entries = l.split("\t")

          // is it an error?
          // yes => is it in a cluster?
          // no => print to file
          entries(12) match {
            case "mis" => {
              numErrors += 1
              
              val piece = entries(2)
              val offset = entries(3).toLong
              val absPos = genome.getAbsPos(piece, offset)

              if (!clusterMembership.contains(absPos)) {
                numErrorsNotInClusters += 1
                
                bw.write(l)
                bw.newLine
              }
            }
            case _ => 
          }
        })
        
	      bw.close

        // output summary
        printf("# errors: %d (%.3f%%)\n", numErrors, (numErrors.toDouble / 100000 * 100.0))
        printf("# errors not in clusters: %d (%.3f%%)\n", numErrorsNotInClusters, (numErrorsNotInClusters.toDouble / numErrors * 100.0))
      }
      case Array("getTruth", alignerSamSubset, trueSam, outFile) => {
        // open output file
        val fw = new java.io.FileWriter(outFile)
        val bw = new java.io.BufferedWriter(fw)
 
        // read in aligner sam subset & record ids contained in the subset
        var source = scala.io.Source.fromFile(alignerSamSubset)
        var lines = source.getLines

        val readIds = Array.fill(100000)(0)
        
        lines.foreach(l => {
          val entries = l.split("\t")
          val idAsInt = masonIdToInt(entries(0))
          readIds(idAsInt) = 1
        })
        
        // read in true sam & write out reads with ids pre-recorded
        source = scala.io.Source.fromFile(trueSam)
        lines = source.getLines
        
        lines.foreach(l => {
          val entries = l.split("\t")
          val idAsInt = masonIdToInt(entries(0))
          if (readIds(idAsInt) == 1) {
            bw.write(l)
            bw.newLine
          }
        })
        
        bw.close
      }
      case Array("explainEasyReadsInClusters", clusterFile, alignerFile, genomeFile, samFile, outFile) => {
        println("Loading genome...")
        val genome = FASTA.read(genomeFile)

        // load cluster file
        // want pos => cluster size
        val clusterSizeByPos = scala.collection.mutable.Map[Long, Int]()
        var source = scala.io.Source.fromFile(clusterFile)
        var lines = source.getLines
        
        println("Loading cluster file...")

        var i = 0
        lines.foreach(l => {
          if (i % 10000000 == 0)
            println("Parsing line " + i + "...")

          // Skip first line
          if (i > 0) {
            val entries = l.split(" ")
            val size = entries(0).toInt
            (1 until entries.length).foreach(n => clusterSizeByPos += ((entries(n).toLong + 1 /* because 0-indexed; want 1-indexed */ , size)))
          }

          i += 1
        })
        
        // load in aligner file
        println("Loading aligner file...")
        source = scala.io.Source.fromFile(alignerFile)
        lines = source.getLines

        val easyReads = Array.fill(100000)(0)

        i = 0
        lines.foreach(l => {
          // skip first line
          if (i > 0) {
            val entries = l.split("\t")
            
            // first column is read id; next 4 columns are codes for the 4 aligners; last column is for cluster membership
            assert(entries.length == 6)
            val numCorrect = entries.toList.slice(1, entries.length - 1).map(_.toInt).count(_ == 1 /* correct */)

            // if in cluster, and all 4 aligners were correct, record id
            if (numCorrect == 4 && entries(entries.length - 1).toInt == 1) {
              easyReads(entries(0).toInt) = 1
              /*
              // works to use snap-based cluster membership determination b/c true pos = aligned pos when all 4 are correct!
              // but, need to read in snap sam so i can get pos (right now, i just have ids)
              clusterSizeByPos.get() match {
                case Some(size) => {
                  // print size to file
                  // update running values so i can compute min, max, average cluster size for "easy" reads
                }
                case None => println("Expected to be in cluster, but isn't.")
              }
              */
            }
          }

          i += 1
        })

	      println("# easy reads: " + easyReads.sum)
        
        // now that i've recorded the ids that correspond to easy reads, traverse the sam file to find out their corresponding cluster sizes
        println("Loading sam file...")
        source = scala.io.Source.fromFile(samFile)
        lines = source.getLines
        
        var minClusterSize = Int.MaxValue
        var maxClusterSize = Int.MinValue
        var runningSum = 0
        var numEasyReads = 0
        
        val fw = new java.io.FileWriter(outFile)
        val bw = new java.io.BufferedWriter(fw)
        
        i = 0
        lines.foreach(l => {
          // if id is recorded as easy/in cluster (1 in corresponding spot in easyReads array), find out its cluster size
          val entries = l.split("\t")
          val idAsInt = masonIdToInt(entries(0))
          if (easyReads(idAsInt) == 1) {
            numEasyReads += 1
            
            val piece = entries(2)
            val offset = entries(3).toLong
            val absPos = genome.getAbsPos(piece, offset)
            
            clusterSizeByPos.get(absPos /* already 1-indexed */) match {
              case Some(size) => {
                if (size < minClusterSize) minClusterSize = size
                if (size > maxClusterSize) maxClusterSize = size
                runningSum += size
                bw.write(size.toString)
                bw.newLine
              }
              case None => println("Alert: expected " + entries(0) + " to be in a cluster.")
            }
          }
          
          i += 1
        })
        
        // print stats
      	println("easyReads.sum: " + easyReads.sum)
      	println("numEasyReads: " + numEasyReads)
        assert(easyReads.sum == numEasyReads)
        println("# easy reads: " + numEasyReads)
        println("Min cluster size: " + minClusterSize)
        println("Max cluster size: " + maxClusterSize)
        println("Average cluster size: " + (runningSum.toDouble / numEasyReads))
        
        bw.close
      }
      case Array("explainHardReadsInClusters", clusterFile, alignerFile, genomeFile, samFile, outFile) => {
        println("Loading genome...")
        val genome = FASTA.read(genomeFile)

        // load cluster file
        // want pos => cluster size
        val clusterSizeByPos = scala.collection.mutable.Map[Long, Int]()
        var source = scala.io.Source.fromFile(clusterFile)
        var lines = source.getLines
        
        println("Loading cluster file...")

        var i = 0
        lines.foreach(l => {
          if (i % 10000000 == 0)
            println("Parsing line " + i + "...")

          // Skip first line
          if (i > 0) {
            val entries = l.split(" ")
            val size = entries(0).toInt
            (1 until entries.length).foreach(n => clusterSizeByPos += ((entries(n).toLong + 1 /* because 0-indexed; want 1-indexed */ , size)))
          }

          i += 1
        })
        
        // load in aligner file
        println("Loading aligner file...")
        source = scala.io.Source.fromFile(alignerFile)
        lines = source.getLines

        val hardReads = Array.fill(100000)(0)

        i = 0
        lines.foreach(l => {
          // skip first line
          if (i > 0) {
            val entries = l.split("\t")
            
            // first column is read id; next 4 columns are codes for the 4 aligners; last column is for cluster membership
            assert(entries.length == 6)
            val numCorrect = entries.toList.slice(1, entries.length - 1).map(_.toInt).count(_ == 1 /* correct */)

            // if in cluster, and at least one aligner did not choose correct answer, record id
            if (numCorrect != 4 && entries(entries.length - 1).toInt == 1)
              hardReads(entries(0).toInt) = 1
          }

          i += 1
        })

	      println("# hard reads: " + hardReads.sum)
        
        // now that i've recorded the ids that correspond to hard reads, traverse the sam file to find out their corresponding cluster sizes
        println("Loading sam file...")
        source = scala.io.Source.fromFile(samFile)
        lines = source.getLines
        
        var minClusterSize = Int.MaxValue
        var maxClusterSize = Int.MinValue
        var runningSum = 0
        var numHardReads = 0
        
        val fw = new java.io.FileWriter(outFile)
        val bw = new java.io.BufferedWriter(fw)
        
        i = 0
        lines.foreach(l => {
          // if id is recorded as hard/in cluster (1 in corresponding spot in hardReads array), find out its cluster size
          val entries = l.split("\t")
          val idAsInt = masonIdToInt(entries(0))
          if (hardReads(idAsInt) == 1) {
            numHardReads += 1
            
            val piece = entries(2)
            val offset = entries(3).toLong
            val absPos = genome.getAbsPos(piece, offset)
            
            clusterSizeByPos.get(absPos /* already 1-indexed */) match {
              case Some(size) => {
                if (size < minClusterSize) minClusterSize = size
                if (size > maxClusterSize) maxClusterSize = size
                runningSum += size
                bw.write(size.toString)
                bw.newLine
              }
              case None => println("Alert: expected " + entries(0) + " to be in a cluster.")
            }
          }
          
          i += 1
        })
        
        // print stats
        assert(hardReads.sum == numHardReads)
        println("# hard reads: " + numHardReads)
        println("Min cluster size: " + minClusterSize)
        println("Max cluster size: " + maxClusterSize)
        println("Average cluster size: " + (runningSum.toDouble / numHardReads))
        
        bw.close
      }
      case Array("prepareAlignerComparisonFile", bowtieFile, bwaFile, novoFile, snapFile, clusterFile, outFile, numReadsStr) => {
        val numReads = numReadsStr.toInt
        val bowtieResults = Array.fill(numReads)(0)
        val bwaResults = Array.fill(numReads)(0)
        val novoResults = Array.fill(numReads)(0)
        val snapResults = Array.fill(numReads)(0)
        val clusterResults = Array.fill(numReads)(0)
        
        // open, parse, & store bowtie results
        var source = scala.io.Source.fromFile(bowtieFile)
        var lines = source.getLines
        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val result = entries(1).toInt
          val i = pairedMasonIdToInt(id)
          bowtieResults(i) = result
        })
        
        // open, parse, & store bwa results
        source = scala.io.Source.fromFile(bwaFile)
        lines = source.getLines
        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val result = entries(1).toInt
          val i = pairedMasonIdToInt(id)
          bwaResults(i) = result
        })

        // open, parse, & store novo results
        source = scala.io.Source.fromFile(novoFile)
        lines = source.getLines
        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val result = entries(1).toInt
          val i = pairedMasonIdToInt(id)
          novoResults(i) = result
        })

        // open, parse, & store snap results
        source = scala.io.Source.fromFile(snapFile)
        lines = source.getLines
        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val result = entries(1).toInt
          val i = pairedMasonIdToInt(id)
          snapResults(i) = result
        })
        
        // open, parse, & store cluster results
        source = scala.io.Source.fromFile(clusterFile)
        lines = source.getLines
        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val result = entries(1).toInt
          val i = pairedMasonIdToInt(id)
          clusterResults(i) = result
        })
        
        // print file
        // format: idAsInt, bowtieResult, bwaResult, novoResult, snapResult, clusterResult
        val fw = new java.io.FileWriter(outFile)
        val bw = new java.io.BufferedWriter(fw)
        
        bw.write(List("idAsInt", "bowtieResult", "bwaResult", "novoResult", "snapResult", "clusterResult").mkString("\t"))
        bw.newLine
        
        (0 until numReads).foreach(r => {
          bw.write(List(r, bowtieResults(r), bwaResults(r), novoResults(r), snapResults(r), clusterResults(r)).mkString("\t"))
          bw.newLine
        })
        
        bw.close
      }
      case Array("prepForClusterLabelingExperiment", genomeFile, clusterFile, alignerFile, samFile) => {
        println("Loading genome...")
        val genome = FASTA.read(genomeFile)

        // load cluster file
        // want pos => cluster id
        var source = scala.io.Source.fromFile(clusterFile)
        var lines = source.getLines
        var i = 0
        val clusterIdByPos = scala.collection.mutable.Map[Long, Long]()
        val clusterSizeById = scala.collection.mutable.Map[Long, Int]()
        
        println("Loading cluster file...")
        lines.foreach(l => {
          if (i % 10000000 == 0)
            println("Parsing line " + i + "...")

          // Skip first line
          if (i > 0) {
            val entries = l.split(" ")
            val size = entries(0).toInt
            val id = entries(1).toLong + 1 // because 0-indexed; want 1-indexed
            clusterSizeById += ((id, size))
            (1 until entries.length).foreach(n => clusterIdByPos += ((entries(n).toLong + 1, id)))  // because 0-indexed; want 1-indexed
          }

          i += 1
        })
        
        // for each read, need to know its aligned pos & how many aligners chose correct pos
        // load in aligner file
        println("Loading aligner file...")
        source = scala.io.Source.fromFile(alignerFile)
        lines = source.getLines

        val readsInClusters = Array.fill(100000 /* this should eventually be passed in */)(-1)
        // -1 => not in cluster
        // 0 => 0 aligners chose correct pos
        // ...
        // 4 => 4 aligners chose correct pos
        
        i = 0
        lines.foreach(l => {
          // skip first line
          if (i > 0) {
            val entries = l.split("\t")
            
            // first column is read id; next 4 columns are codes for the 4 aligners; last column is for cluster membership
            assert(entries.length == 6)
            val numCorrect = entries.toList.slice(1, entries.length - 1).map(_.toInt).count(_ == 1 /* correct */)

            // if in cluster, record id and how many aligners were correct
            if (entries(entries.length - 1).toInt == 1)
              readsInClusters(entries(0).toInt) = numCorrect
          }

          i += 1
        })


        // traverse sam file & if a read is in a cluster (use readsInClusters), use clusterIdByPos to figure out which
        val readsHittingClustersById = scala.collection.mutable.Map[Long, Array[Int]]()
        
        println("Loading sam file...")
        source = scala.io.Source.fromFile(samFile)
        lines = source.getLines
        
        i = 0
        
        lines.foreach(l => {
          val entries = l.split("\t")
          val idAsInt = masonIdToInt(entries(0))
          val readDifficultyScore = readsInClusters(idAsInt)
          if (readDifficultyScore > -1) {
            val piece = entries(2)
            val offset = entries(3).toLong
            val absPos = genome.getAbsPos(piece, offset)
            
            clusterIdByPos.get(absPos /* already 1-indexed */) match {
              case Some(id) => {
                readsHittingClustersById.get(id) match {
                  case Some(readsArray) => {
                    readsArray(readDifficultyScore) += 1
                    readsHittingClustersById.update(id, readsArray)
                  }
                  case None => {
                    val readsArray = Array.fill(5 /* 0,1,2,3,4 */)(0)
                    readsArray(readDifficultyScore) += 1
                    readsHittingClustersById += ((id, readsArray))
                  }
                }
              }
              case None => println("Alert: expected " + entries(0) + " to be in a cluster.")
            }
          }
          
          i += 1
        })
        
        println("# clusters with more than one read: " + readsHittingClustersById.values.count(_.sum > 1))
        println("# mixed clusters: " + readsHittingClustersById.values.filter(i => i.count(_ != 0) > 1).toList.length)
        
        // for each of the mixed clusters, get size & # of each category in that cluster (eg 0-5, 1-2, ..., 4-1)
        println("clusterId, clusterSize, 0, 1, 2, 3, 4")
        //res0.keys.filter(k => res0.get(k).get.sum > 1).foreach(k => println(res0.get(k).get.mkString(",")))
        //println(List(1,2).mkString(",") + "," + Array(3,4,5).mkString(","))
        readsHittingClustersById.keys.filter(k => readsHittingClustersById.get(k).get.count(_ != 0) > 1)
          .foreach(k => println(List(k, clusterSizeById.get(k).get).mkString(", ") + ", " + readsHittingClustersById.get(k).get.mkString(",")))
      }
      case Array("createClusterIdByPosMap", clusterFile, outFile) => {
        // load cluster file
        // want pos => cluster id
        val clusterIdByPos = scala.collection.mutable.Map[Long, Long]() // want to figure out how to read/write this to file
        var source = scala.io.Source.fromFile(clusterFile)
        var lines = source.getLines
        var i = 0
        
        println("Loading cluster file...")
        lines.foreach(l => {
          if (i % 10000000 == 0)
            println("Parsing line " + i + "...")

          // Skip first line
          if (i > 0) {
            val entries = l.split(" ")
            //val size = entries(0).toInt
            val id = entries(1).toLong + 1 /* because 0-indexed; want 1-indexed */
            (1 until entries.length).foreach(n => clusterIdByPos += ((entries(n).toLong + 1 /* because 0-indexed; want 1-indexed */, id)))
          }

          i += 1
        })

        val fos = new java.io.FileOutputStream(outFile)
        val oos = new java.io.ObjectOutputStream(fos)
        oos.writeObject(clusterIdByPos)
        oos.close
      }
      case Array("correlateAlignerAgreementWithCorrectness", bowtieSam, bwaSam, novoSam, snapSam, genomeFile, numReadsStr, toleranceStr) => {
        // for each read, check:  do all 4 aligners agree?  if so, are they correct?

        println("Loading genome...")
        val genome = FASTA.read(genomeFile)
        println("Done loading genome.")
        
        // for each read, store ID, pos chosen by aligner, whether it was correct
        val numReads = numReadsStr.toInt
        val bowtiePos = Array.fill(numReads)(0L)
        val bowtieResult = Array.fill(numReads)(0)
        val bwaPos = Array.fill(numReads)(0L)
        val bwaResult = Array.fill(numReads)(0)
        val novoPos = Array.fill(numReads)(0L)
        val novoResult = Array.fill(numReads)(0)
        val snapPos = Array.fill(numReads)(0L)
        val snapResult = Array.fill(numReads)(0)
        
        // open, parse, & store bowtie results
        println("Parsing Bowtie2 file...")
        var source = scala.io.Source.fromFile(bowtieSam)
        var lines = source.getLines
        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val i = pairedMasonIdToInt(id)

          val piece = entries(2)
          val offset = entries(3).toLong
          val absPos = genome.getAbsPos(piece, offset)

          // check to see if it was aligned correctly
          val code = 
          entries(12) match {
            case "mis" => -1
            case "un" => 0
            case _ => 1
          }

          bowtiePos(i) = absPos
          bowtieResult(i) = code
        })

        // open, parse, & store bwa results
        println("Parsing BWA file...")
        source = scala.io.Source.fromFile(bwaSam)
        lines = source.getLines
        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val i = pairedMasonIdToInt(id)

          val piece = entries(2)
          val offset = entries(3).toLong
          val absPos = genome.getAbsPos(piece, offset)

          // check to see if it was aligned correctly
          val code = 
          entries(12) match {
            case "mis" => -1
            case "un" => 0
            case _ => 1
          }

          bwaPos(i) = absPos
          bwaResult(i) = code
        })
        
        // open, parse, & store novoalign results
        println("Parsing Novoalign file...")
        source = scala.io.Source.fromFile(novoSam)
        lines = source.getLines
        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val i = pairedMasonIdToInt(id)

          val piece = entries(2)
          val offset = entries(3).toLong
          val absPos = genome.getAbsPos(piece, offset)

          // check to see if it was aligned correctly
          val code = 
          entries(12) match {
            case "mis" => -1
            case "un" => 0
            case _ => 1
          }

          novoPos(i) = absPos
          novoResult(i) = code
        })
        
        // open, parse, & store snap results
        println("Parsing SNAP file...")
        source = scala.io.Source.fromFile(snapSam)
        lines = source.getLines
        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val i = pairedMasonIdToInt(id)

          val piece = entries(2)
          val offset = entries(3).toLong
          val absPos = genome.getAbsPos(piece, offset)

          // check to see if it was aligned correctly
          val code = 
          entries(12) match {
            case "mis" => -1
            case "un" => 0
            case _ => 1
          }

          snapPos(i) = absPos
          snapResult(i) = code
        })
        
        // check agreement
        println("Checking agreement...")
        var numReadsAgree = 0
        var numReadsAgreeCorrect = 0
        val tolerance = toleranceStr.toInt
        
        (0 until numReads).foreach(i => {
          // do aligners agree?
          val minPos = List(bowtiePos(i), bwaPos(i), novoPos(i), snapPos(i)).min
          val maxPos = List(bowtiePos(i), bwaPos(i), novoPos(i), snapPos(i)).max
          if ((maxPos - minPos) < tolerance) {
            numReadsAgree += 1

            // are all 4 aligners correct?
            if (bowtieResult(i) == 1 && bwaResult(i) == 1 && novoResult(i) == 1 && snapResult(i) == 1)
              numReadsAgreeCorrect += 1
          }
        })
        
        printf("Aligners agreed on %d reads.  Out of those, %d (%.3f%%) were correct.\n", numReadsAgree, numReadsAgreeCorrect, 
          (numReadsAgreeCorrect.toDouble / numReadsAgree * 100.0))
      }
      case Array("checkRealDataAlignerAgreementWrtClusters", bowtieSam, bwaSam, novoSam, snapSam, clusterFile, genomeFile, /*numReadsStr*/numPairsStr, toleranceStr) => {
        println("Loading genome...")
        val genome = FASTA.read(genomeFile)
        println("Done loading genome.")

        println("Loading clusters...")
        val clusterMembership = getClusterMembers(clusterFile)
        println("Done loading clusters.")
        
        /*
        val numReads = numReadsStr.toInt
        val bowtiePos = Array.fill(numReads)(0L)
        val bwaPos = Array.fill(numReads)(0L)
        val novoPos = Array.fill(numReads)(0L)
        val snapPos = Array.fill(numReads)(0L)
        */
        val numPairs = numPairsStr.toInt // assumes < 2B pairs for whole genome
        val bowtiePos = Array.fill(2)(Array.fill(numPairs)(0L))
        val bwaPos = Array.fill(2)(Array.fill(numPairs)(0L))
        val novoPos = Array.fill(2)(Array.fill(numPairs)(0L))
        val snapPos = Array.fill(2)(Array.fill(numPairs)(0L))

        // open, parse, & store bowtie results
        println("Parsing Bowtie2 file...")
        var samEntries = SAM.read(bowtieSam)
        var i = 0
        samEntries.foreach(samEntry => {
          // get id
          // assumes id is in format like ERR091571.1
          val idEntries = samEntry.readId.split("\\.")
          val readNum = idEntries(1).toInt
          //val index = 2 * readNum + (i % 2)

          val absPos = genome.getAbsPos(samEntry.piece, samEntry.position)
          // infer if it's first or 2nd based on whether i is even or odd
          // i is even => first read; i is odd => second read
          bowtiePos(i % 2)(readNum - 1) = absPos
          
          i += 1
        })

        // open, parse, & store bwa results
        println("Parsing BWA file...")
        samEntries = SAM.read(bwaSam)
        i = 0
        samEntries.foreach(samEntry => {
          // get id
          // assumes id is in format like ERR091571.1
          val idEntries = samEntry.readId.split("\\.")
          val readNum = idEntries(1).toInt
          //val index = 2 * readNum + (i % 2)

          val absPos = genome.getAbsPos(samEntry.piece, samEntry.position)
          // infer if it's first or 2nd based on whether i is even or odd
          // i is even => first read; i is odd => second read
          bwaPos(i % 2)(readNum - 1) = absPos
          
          i += 1
        })

        // open, parse, & store novo results
        println("Parsing Novoalign file...")
        samEntries = SAM.read(novoSam)
        i = 0
        samEntries.foreach(samEntry => {
          // get id
          // assumes id is in format like ERR091571.1
          val idEntries = samEntry.readId.split("\\.")
          val readNum = idEntries(1).toInt
          //val index = 2 * readNum + (i % 2)

          val absPos = genome.getAbsPos(samEntry.piece, samEntry.position)
          // infer if it's first or 2nd based on whether i is even or odd
          // i is even => first read; i is odd => second read
          novoPos(i % 2)(readNum - 1) = absPos
          
          i += 1
        })
        
        // open, parse, & store snap results
        println("Parsing SNAP file...")
        samEntries = SAM.read(snapSam)
        i = 0
        samEntries.foreach(samEntry => {
          // get id
          // assumes id is in format like ERR091571.1
          val idEntries = samEntry.readId.split("\\.")
          val readNum = idEntries(1).toInt
          //val index = 2 * readNum + (i % 2)

          val absPos = genome.getAbsPos(samEntry.piece, samEntry.position)
          // infer if it's first or 2nd based on whether i is even or odd
          // i is even => first read; i is odd => second read
          snapPos(i % 2)(readNum - 1) = absPos
          
          i += 1
        })
        
        // check agreement
        println("Checking agreement...")
        var numReadsAgree = 0L
        var numReadsAgreeInClusters = 0L
        val tolerance = toleranceStr.toInt
        
        (0 until numPairs).foreach(i => {
          List(0, 1).foreach(whichRead => {
            // do aligners agree?
            val minPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).min
            val maxPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).max
            if ((maxPos - minPos) < tolerance) {
              numReadsAgree += 1 

              // is read in cluster?  count YES if any of the aligners' chosen positions are in clusters
              if (clusterMembership.contains(bowtiePos(whichRead)(i)) || 
                clusterMembership.contains(bwaPos(whichRead)(i)) || 
                clusterMembership.contains(novoPos(whichRead)(i)) || 
                clusterMembership.contains(snapPos(whichRead)(i)))
                numReadsAgreeInClusters += 1
            }
          })
        })
        
        printf("Aligners agreed on %d reads.  Out of those, %d (%.3f%%) were in clusters.\n", numReadsAgree, numReadsAgreeInClusters, 
          (numReadsAgreeInClusters.toDouble / numReadsAgree * 100.0))
        
        // see how many aligners agree
        /*
        val numAlignersAgree = Array.fill(numReads)(0)
        val inCluster = Array.fill(numReads)(0)
        */
        val numAlignersAgree = Array.fill(2)(Array.fill(numPairs)(0))
        val inCluster = Array.fill(2)(Array.fill(numPairs)(0))
        var done = false
        
        (0 until numPairs).foreach(i => {
          done = false
          List(0, 1).foreach(whichRead => {
            // do all 4 aligners agree?
            var minPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).min
            var maxPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).max
            if ((maxPos - minPos) < tolerance) {
              numAlignersAgree(whichRead)(i) = 4
              done = true
            }

            // do 3 aligners agree?
            // check bowtie, bwa, novo
            minPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), novoPos(whichRead)(i)).min
            maxPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), novoPos(whichRead)(i)).max
            if (!done && (maxPos - minPos) < tolerance) {
              numAlignersAgree(whichRead)(i) = 3
              done = true
            }

            // check bowtie, bwa, snap
            minPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), snapPos(whichRead)(i)).min
            maxPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), snapPos(whichRead)(i)).max
            if (!done && (maxPos - minPos) < tolerance) {
              numAlignersAgree(whichRead)(i) = 3
              done = true
            }

            // check bowtie, novo, snap
            minPos = List(bowtiePos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).min
            maxPos = List(bowtiePos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).max
            if (!done && (maxPos - minPos) < tolerance) {
              numAlignersAgree(whichRead)(i) = 3
              done = true
            }

            // check bwa, novo, snap
            minPos = List(bwaPos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).min
            maxPos = List(bwaPos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).max
            if (!done && (maxPos - minPos) < tolerance) {
              numAlignersAgree(whichRead)(i) = 3
              done = true
            }

            // do 2 aligners agree?
            // check bowtie, bwa
            if (!done && math.abs(bowtiePos(whichRead)(i) - bwaPos(whichRead)(i)) < tolerance) {
              numAlignersAgree(whichRead)(i) = 2
              done = true
            }

            // check bowtie, novo
            if (!done && math.abs(bowtiePos(whichRead)(i) - novoPos(whichRead)(i)) < tolerance) {
              numAlignersAgree(whichRead)(i) = 2
              done = true
            }

            // check bowtie, snap
            if (!done && math.abs(bowtiePos(whichRead)(i) - snapPos(whichRead)(i)) < tolerance) {
              numAlignersAgree(whichRead)(i) = 2
              done = true
            }

            // check bwa, novo
            if (!done && math.abs(bwaPos(whichRead)(i) - novoPos(whichRead)(i)) < tolerance) {
              numAlignersAgree(whichRead)(i) = 2
              done = true
            }

            // check bwa, snap
            if (!done && math.abs(bwaPos(whichRead)(i) - snapPos(whichRead)(i)) < tolerance) {
              numAlignersAgree(whichRead)(i) = 2
              done = true
            }

            // check novo, snap
            if (!done && math.abs(novoPos(whichRead)(i) - snapPos(whichRead)(i)) < tolerance) {
              numAlignersAgree(whichRead)(i) = 2
              done = true
            }

            done = false

            // check if it's in a cluster
            // def:  a read is in a cluster if ANY of the aligners place it there
            if (clusterMembership.contains(bowtiePos(whichRead)(i)) || 
              clusterMembership.contains(bwaPos(whichRead)(i)) || 
              clusterMembership.contains(novoPos(whichRead)(i)) || 
              clusterMembership.contains(snapPos(whichRead)(i)))
              inCluster(whichRead)(i) = 1          
          })
        })

        val inClustersNAgree = Array.fill(5)(0L)
        val notInClustersNAgree = Array.fill(5)(0L)
        
        (0 until numPairs).foreach(i => {
          List(0, 1).foreach(whichRead => {
            val numAgree = numAlignersAgree(whichRead)(i)
            if (inCluster(whichRead)(i) == 1)
              inClustersNAgree(numAgree) += 1
            else
              notInClustersNAgree(numAgree) += 1
          })
        })
        
        // print agreement table
        println("#Aligners\tInClusters\tNotInClusters\tPercentInClusters")
        // skip one b/c it doesn't make sense
        List(0, 2, 3, 4).foreach(i => {
          printf("%d\t%d\t%d\t%.3f%%\n", i, inClustersNAgree(i), notInClustersNAgree(i), (inClustersNAgree(i).toDouble / (inClustersNAgree(i) + notInClustersNAgree(i)) * 100.0))
        })
      }
      case Array("similarRegionsToBedFormat", genomeFile, clusterFile, outFile, readLenStr) => {
        val readLen = readLenStr.toInt
        
        println("Loading genome...")
        val genome = FASTA.read(genomeFile)
        println("Done loading genome.")

        // Open cluster file
        val source = scala.io.Source.fromFile(clusterFile)
        val lines = source.getLines
        
        // Open output file
        val fw = new java.io.FileWriter(outFile)
        val bw = new java.io.BufferedWriter(fw)
                
        // Parse cluster file.
        // For each position in a cluster (ie, for each similar region), output the corresponding bed format entry.
        var i = 0
        lines.foreach(l => {
          if (i % 10000000 == 0)
            println("Parsing line " + i + "...")

          // Skip first line
          if (i > 0) {
            val entries = l.split(" ")
            // entries(0) is cluster size
            (1 until entries.length).foreach(n => {
              val absPos = entries(n).toLong // 0-indexed
              val (piece, pos) = genome.getLocation(absPos) // 0-indexed
              // format:  chrName \t startPos \t endPos
              // [startPos, endPos) -- startPos is included in the region, while endPos is not
              // thus, to get start & end pos:  start pos is 0-indexed start pos (obtained from getLocation); endPos = startPos + readLen
              bw.write(piece + "\t" + pos.toString + "\t" + (pos + readLen).toString)
              bw.newLine
            })
          }

          i += 1
        })
        
        bw.close        
      }
      case Array("analyzeContiguousSimilarRegions", clusterFile, genomeFile, readLenStr, similarRegionLengthFile, interSimilarRegionLengthFile) => {
        val readLen = readLenStr.toInt
        
        // load genome
        println("Loading genome...")
        val genome = FASTA.read(genomeFile)
        println("Done loading genome.")
        
        // create bitmap for genome
        // set any pos to true that is within a substring that is in a similar region
        val genomeSize = genome.totalSize
        val midPoint = (genomeSize / 2).toInt

        // 0 => not in similar region
        // 1 => in similar region
        // keep everything 0-indexed (both genome and similar regions are)
        val bitmapL = Array.fill(midPoint)(0) // 0 to midPoint - 1
        val bitmapH = Array.fill(midPoint)(0) // midPoint to genomeSize - 1
        
        // for each substring in similar regions, set bits to 1
        val source = scala.io.Source.fromFile(clusterFile)
        val lines = source.getLines
        var j = 0
        lines.foreach(l => {
          if (j % 10000000 == 0)
            println("Parsing line " + j + "...")
          
            // Skip first line
            if (j > 0) {
              val entries = l.split(" ")
              // entries(0) is cluster size
              (1 until entries.length).foreach(n => {
                val startPos = entries(n).toLong
                val endPos = startPos + readLen - 1
                
                if (startPos < midPoint && endPos < midPoint) {
                  (startPos to endPos).foreach(p => bitmapL(p.toInt) = 1)
                } else if (startPos < midPoint && endPos >= midPoint) {
                  (startPos until midPoint).foreach(p => bitmapL(p.toInt) = 1)
                  (midPoint.toLong to endPos).foreach(p => bitmapH((p - midPoint).toInt) = 1)
                } else if (startPos >= midPoint && endPos >= midPoint) {
                  (startPos to endPos).foreach(p => bitmapH((p - midPoint).toInt) = 1)
                }
              })
            }

            j += 1
        })
        
        // analyze bitmap
        // Q:  how many bits are set to 1?
        val numPosInSimilarRegions = bitmapL.sum + bitmapH.sum
        printf("Pos in similar regions: %d/%d (%.3f%%)\n", numPosInSimilarRegions, genomeSize, (numPosInSimilarRegions.toDouble / genomeSize * 100.0))
        
        // Q:  how many distinct similar regions are there?
        // Q:  what is their length distribution?
        // Q:  what is their inter-length distribution?
        var numSimilarRegions = 0
        var prevPos = 0
        var similarRegionLengths: List[Long] = Nil
        var interSimilarRegionLengths: List[Long] = Nil
        var startPosOfLastSimilarRegion = 0L
        var endPosOfLastSimilarRegion = 0L
        var i = 0L
        while (i < genomeSize) {
          val bitmapVal = 
            if (i < midPoint) bitmapL(i.toInt)
            else bitmapH((i - midPoint).toInt)

          // if prevPos == bitmapVal == 0, do nothing
          if (prevPos == 0 && (bitmapVal == 1 || i == (genomeSize - 1))) {
            // you've reached the beginning of a similar region
            numSimilarRegions += 1
            prevPos = 1
            startPosOfLastSimilarRegion = i
            interSimilarRegionLengths = 
	            if (numSimilarRegions == 1 || i == (genomeSize - 1))	// boundary case
	              if (i - endPosOfLastSimilarRegion > 0)
	    	          (i - endPosOfLastSimilarRegion) :: interSimilarRegionLengths
		            else
		              interSimilarRegionLengths
		          else
		            (i - endPosOfLastSimilarRegion - 1) :: interSimilarRegionLengths
          }
          else if (prevPos == 1 && (bitmapVal == 0 || i == (genomeSize - 1))) {
            // you've reached the end of a similar region
            prevPos = 0
            similarRegionLengths = 
	            if (i == (genomeSize - 1))
	              (i - startPosOfLastSimilarRegion + 1) :: similarRegionLengths
	            else
	              (i - startPosOfLastSimilarRegion) :: similarRegionLengths
            endPosOfLastSimilarRegion = i - 1
          }
          // if prevPos == bitmapVal == 1, do nothing
          
          i += 1
        }

        println("# distinct similar regions: " + numSimilarRegions)
        // output similar region lengths list to file
        var fw = new java.io.FileWriter(similarRegionLengthFile)
        var bw = new java.io.BufferedWriter(fw)
        
        similarRegionLengths.foreach(l => {
          bw.write(l.toString)
          bw.newLine
        })
        
        bw.close
        
        // output inter similar region lengths list to file
        fw = new java.io.FileWriter(interSimilarRegionLengthFile)
        bw = new java.io.BufferedWriter(fw)
        
        interSimilarRegionLengths.foreach(l => {
          bw.write(l.toString)
          bw.newLine
        })
        
        bw.close
      }
      case Array("checkSimulatedDataAlignerAgreementWrtClusters", bowtieSam, bwaSam, novoSam, snapSam, clusterFile, genomeFile, numPairsStr, toleranceStr) => {
        println("Loading genome...")
        val genome = FASTA.read(genomeFile)
        println("Done loading genome.")

        println("Loading clusters...")
        val clusterMembership = getClusterMembers(clusterFile)
        println("Done loading clusters.")
        
        val numPairs = numPairsStr.toInt // assumes < 2B pairs for whole genome
        val bowtiePos = Array.fill(2)(Array.fill(numPairs)(0L))
        val bwaPos = Array.fill(2)(Array.fill(numPairs)(0L))
        val novoPos = Array.fill(2)(Array.fill(numPairs)(0L))
        val snapPos = Array.fill(2)(Array.fill(numPairs)(0L))

        /*
        var source = scala.io.Source.fromFile(bowtieFile)
        var lines = source.getLines
        lines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0)
          val result = entries(1).toInt
          val i = pairedMasonIdToInt(id)
          bowtieResults(i) = result
        })
        */

        // open, parse, & store bowtie results
        println("Parsing Bowtie2 file...")
        var source = scala.io.Source.fromFile(bowtieSam)
        var lines = source.getLines
        var i = 0
        lines.foreach(l => {
          // get id
          val entries = l.split("\t")
          val id = entries(0)
          val idEntries = id.split("/")
          val whichRead = idEntries(1).toInt // 1 or 2
          val idPrefixEntries = idEntries(0).split("\\.")
          val readNum = idPrefixEntries(idPrefixEntries.length - 1).toInt // [0, numPairs)
          val absPos = genome.getAbsPos(entries(2), entries(3).toLong)
          bowtiePos(whichRead - 1)(readNum) = absPos
          
          i += 1
        })

        // open, parse, & store bwa results
        println("Parsing BWA file...")
        source = scala.io.Source.fromFile(bwaSam)
        lines = source.getLines
        i = 0
        lines.foreach(l => {
          // get id
          val entries = l.split("\t")
          val id = entries(0)
          val idEntries = id.split("/")
          val whichRead = idEntries(1).toInt // 1 or 2
          val idPrefixEntries = idEntries(0).split("\\.")
          val readNum = idPrefixEntries(idPrefixEntries.length - 1).toInt // [0, numPairs)
          val absPos = genome.getAbsPos(entries(2), entries(3).toLong)
          bwaPos(whichRead - 1)(readNum) = absPos
          
          i += 1
        })

        // open, parse, & store novo results
        println("Parsing Novoalign file...")
        source = scala.io.Source.fromFile(novoSam)
        lines = source.getLines
        i = 0
        lines.foreach(l => {
          // get id
          val entries = l.split("\t")
          val id = entries(0)
          val idEntries = id.split("/")
          val whichRead = idEntries(1).toInt // 1 or 2
          val idPrefixEntries = idEntries(0).split("\\.")
          val readNum = idPrefixEntries(idPrefixEntries.length - 1).toInt // [0, numPairs)
          val absPos = genome.getAbsPos(entries(2), entries(3).toLong)
          novoPos(whichRead - 1)(readNum) = absPos
          
          i += 1
        })
        
        // open, parse, & store snap results
        println("Parsing SNAP file...")
        source = scala.io.Source.fromFile(snapSam)
        lines = source.getLines
        i = 0
        lines.foreach(l => {
          // get id
          val entries = l.split("\t")
          val id = entries(0)
          val idEntries = id.split("/")
          val whichRead = idEntries(1).toInt // 1 or 2
          val idPrefixEntries = idEntries(0).split("\\.")
          val readNum = idPrefixEntries(idPrefixEntries.length - 1).toInt // [0, numPairs)
          val absPos = genome.getAbsPos(entries(2), entries(3).toLong)
          snapPos(whichRead - 1)(readNum) = absPos
          
          i += 1
        })
        
        // check agreement
        println("Checking agreement...")
        var numReadsAgree = 0L
        var numReadsAgreeInClusters = 0L
        val tolerance = toleranceStr.toInt
        
        (0 until numPairs).foreach(i => {
          List(0, 1).foreach(whichRead => {
            // do aligners agree?
            val minPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).min
            val maxPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).max
            if ((maxPos - minPos) < tolerance) {
              numReadsAgree += 1 

              // is read in cluster?  count YES if any of the aligners' chosen positions are in clusters
              if (clusterMembership.contains(bowtiePos(whichRead)(i)) || 
                clusterMembership.contains(bwaPos(whichRead)(i)) || 
                clusterMembership.contains(novoPos(whichRead)(i)) || 
                clusterMembership.contains(snapPos(whichRead)(i)))
                numReadsAgreeInClusters += 1
            }
          })
        })
        
        printf("Aligners agreed on %d reads.  Out of those, %d (%.3f%%) were in clusters.\n", numReadsAgree, numReadsAgreeInClusters, 
          (numReadsAgreeInClusters.toDouble / numReadsAgree * 100.0))
        
        // see how many aligners agree
        /*
        val numAlignersAgree = Array.fill(numReads)(0)
        val inCluster = Array.fill(numReads)(0)
        */
        val numAlignersAgree = Array.fill(2)(Array.fill(numPairs)(0))
        val inCluster = Array.fill(2)(Array.fill(numPairs)(0))
        var done = false
        
        (0 until numPairs).foreach(i => {
          done = false
          List(0, 1).foreach(whichRead => {
            // do all 4 aligners agree?
            var minPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).min
            var maxPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).max
            if ((maxPos - minPos) < tolerance) {
              numAlignersAgree(whichRead)(i) = 4
              done = true
            }

            // do 3 aligners agree?
            // check bowtie, bwa, novo
            minPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), novoPos(whichRead)(i)).min
            maxPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), novoPos(whichRead)(i)).max
            if (!done && (maxPos - minPos) < tolerance) {
              numAlignersAgree(whichRead)(i) = 3
              done = true
            }

            // check bowtie, bwa, snap
            minPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), snapPos(whichRead)(i)).min
            maxPos = List(bowtiePos(whichRead)(i), bwaPos(whichRead)(i), snapPos(whichRead)(i)).max
            if (!done && (maxPos - minPos) < tolerance) {
              numAlignersAgree(whichRead)(i) = 3
              done = true
            }

            // check bowtie, novo, snap
            minPos = List(bowtiePos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).min
            maxPos = List(bowtiePos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).max
            if (!done && (maxPos - minPos) < tolerance) {
              numAlignersAgree(whichRead)(i) = 3
              done = true
            }

            // check bwa, novo, snap
            minPos = List(bwaPos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).min
            maxPos = List(bwaPos(whichRead)(i), novoPos(whichRead)(i), snapPos(whichRead)(i)).max
            if (!done && (maxPos - minPos) < tolerance) {
              numAlignersAgree(whichRead)(i) = 3
              done = true
            }

            // do 2 aligners agree?
            // check bowtie, bwa
            if (!done && math.abs(bowtiePos(whichRead)(i) - bwaPos(whichRead)(i)) < tolerance) {
              numAlignersAgree(whichRead)(i) = 2
              done = true
            }

            // check bowtie, novo
            if (!done && math.abs(bowtiePos(whichRead)(i) - novoPos(whichRead)(i)) < tolerance) {
              numAlignersAgree(whichRead)(i) = 2
              done = true
            }

            // check bowtie, snap
            if (!done && math.abs(bowtiePos(whichRead)(i) - snapPos(whichRead)(i)) < tolerance) {
              numAlignersAgree(whichRead)(i) = 2
              done = true
            }

            // check bwa, novo
            if (!done && math.abs(bwaPos(whichRead)(i) - novoPos(whichRead)(i)) < tolerance) {
              numAlignersAgree(whichRead)(i) = 2
              done = true
            }

            // check bwa, snap
            if (!done && math.abs(bwaPos(whichRead)(i) - snapPos(whichRead)(i)) < tolerance) {
              numAlignersAgree(whichRead)(i) = 2
              done = true
            }

            // check novo, snap
            if (!done && math.abs(novoPos(whichRead)(i) - snapPos(whichRead)(i)) < tolerance) {
              numAlignersAgree(whichRead)(i) = 2
              done = true
            }

            done = false

            // check if it's in a cluster
            // def:  a read is in a cluster if ANY of the aligners place it there
            if (clusterMembership.contains(bowtiePos(whichRead)(i)) || 
              clusterMembership.contains(bwaPos(whichRead)(i)) || 
              clusterMembership.contains(novoPos(whichRead)(i)) || 
              clusterMembership.contains(snapPos(whichRead)(i)))
              inCluster(whichRead)(i) = 1          
          })
        })

        val inClustersNAgree = Array.fill(5)(0L)
        val notInClustersNAgree = Array.fill(5)(0L)
        
        (0 until numPairs).foreach(i => {
          List(0, 1).foreach(whichRead => {
            val numAgree = numAlignersAgree(whichRead)(i)
            if (inCluster(whichRead)(i) == 1)
              inClustersNAgree(numAgree) += 1
            else
              notInClustersNAgree(numAgree) += 1
          })
        })
        
        // print agreement table
        println("#Aligners\tInClusters\tNotInClusters\tPercentInClusters")
        // skip one b/c it doesn't make sense
        List(0, 2, 3, 4).foreach(i => {
          printf("%d\t%d\t%d\t%.3f%%\n", i, inClustersNAgree(i), notInClustersNAgree(i), (inClustersNAgree(i).toDouble / (inClustersNAgree(i) + notInClustersNAgree(i)) * 100.0))
        })
      }
      case Array("computeSimilarRegionsBitmap", clusterFile, genomeFile, readLenStr, bitmapOutFile) => {
        val readLen = readLenStr.toInt
        
        // load genome
        println("Loading genome...")
        val genome = FASTA.read(genomeFile)
        println("Done loading genome.")
        
        // create bitmap for genome
        // set any pos to true that is within a substring that is in a similar region
        val genomeSize = genome.totalSize
        val midPoint = (genomeSize / 2).toInt

        // 0 => not in similar region
        // 1 => in similar region
        // keep everything 0-indexed (both genome and similar regions are)
        val bitmapL = Array.fill(midPoint)(0) // 0 to midPoint - 1
        val bitmapH = Array.fill(midPoint)(0) // midPoint to genomeSize - 1
        
        // for each substring in similar regions, set bits to 1
        val source = scala.io.Source.fromFile(clusterFile)
        val lines = source.getLines
        var j = 0
        lines.foreach(l => {
          if (j % 10000000 == 0)
            println("Parsing line " + j + "...")
          
            // Skip first line
            if (j > 0) {
              val entries = l.split(" ")
              // entries(0) is cluster size
              (1 until entries.length).foreach(n => {
                val startPos = entries(n).toLong
                val endPos = startPos + readLen - 1
                
                if (startPos < midPoint && endPos < midPoint) {
                  (startPos to endPos).foreach(p => bitmapL(p.toInt) = 1)
                } else if (startPos < midPoint && endPos >= midPoint) {
                  (startPos until midPoint).foreach(p => bitmapL(p.toInt) = 1)
                  (midPoint.toLong to endPos).foreach(p => bitmapH((p - midPoint).toInt) = 1)
                } else if (startPos >= midPoint && endPos >= midPoint) {
                  (startPos to endPos).foreach(p => bitmapH((p - midPoint).toInt) = 1)
                }
              })
            }

            j += 1
        })
        
        // Save bitmap to file
        println("Saving bitmap to file...")
        val fouts = new FileOutputStream(bitmapOutFile);

        // Write object with ObjectOutputStream
        val outobj= new ObjectOutputStream (fouts);

        // Write object out to disk
        outobj.writeObject ( bitmapL );
        outobj.writeObject ( bitmapH );
        
        outobj.close
        
        /*
        // Try reading it back in
        println("Reading bitmap from file...")
        val infile = new FileInputStream(bitmapOutFile);

        val inobj = new ObjectInputStream (infile);

        val objL = inobj.readObject();
        val objH = inobj.readObject();

        val bitmapL2 = objL.asInstanceOf[Array[Int]]
        val bitmapH2 = objH.asInstanceOf[Array[Int]]
        
        println("bitmapL sum: " + bitmapL2.sum)
        println("bitmapH sum: " + bitmapH2.sum)
        */
      }
      case _ => println("Incorrect args.")
    }
  }

  def getClusterMembers(clusterFilename: String): scala.collection.mutable.Set[Long] = {
    val source = scala.io.Source.fromFile(clusterFilename)
    val lines = source.getLines
    val clusterMembership = scala.collection.mutable.Set[Long]()

    var i = 0
    lines.foreach(l => {
      if (i % 10000000 == 0)
        println("Parsing line " + i + "...")

      // Skip first line
      if (i > 0) {
        val entries = l.split(" ")
        // entries(0) is cluster size
        (1 until entries.length).foreach(n => clusterMembership += (entries(n).toLong + 1 /* because 0-indexed; want 1-indexed */ ))
      }
        
      i += 1
    })

    clusterMembership
  }

  def getClusterIdByPos(clusterFile: String, want0Indexed: Boolean) = {
    val source = scala.io.Source.fromFile(clusterFile)
    val lines = source.getLines
    var i = 0L
    val clusterIdByPos = scala.collection.mutable.Map[Long, Long]()
    
    println("Loading cluster file...")
    lines.foreach(l => {
      if (i % 10000000 == 0)
        println("Parsing line " + i + "...")

      // Skip first line
      if (i > 0) {
        val entries = l.split(" ")
        val size = entries(0).toInt
        //val id = entries(1).toLong + 1 // because 0-indexed; want 1-indexed
        //(1 until entries.length).foreach(n => clusterIdByPos += ((entries(n).toLong + 1, id)))  // because 0-indexed; want 1-indexed
        val id = 
          if (want0Indexed)
            entries(1).toLong // 0-indexed
          else
            entries(1).toLong + 1 // 1-indexed
        (1 until entries.length).foreach(n => {
          val pos = 
          if (want0Indexed)
            entries(n).toLong // 0-indexed
          else
            entries(n).toLong + 1 // 1-indexed
          clusterIdByPos += ((pos, id))
        })
      }

      i += 1
    })
    
    clusterIdByPos
  }

  def masonIdToInt(id: String): Int = {
    // assumes id formatted like /scratch/mason/hg19.illumina.hn2.sq.N100000.n100.fq.000000008/1
    // will pull out # at the end -- here, 8
    
    val entries = id.split("/")
    val entries2 = entries(3).split("\\.")
    entries2(entries2.length - 1).toInt
  }
  
  def pairedMasonIdToInt(id: String): Int = {
    // assumes id formatted like hg19.illumina.hn2.sq.mp.ll375.le100.N5000000.n100.s0.fq.000000000/1 or /2
    // will pull out # at the end -- here 0 -- and which read in pair this is -- here 1
    
    val entries = id.split("/")
    val entries2 = entries(entries.length - 2).split("\\.")
    2 * entries2(entries2.length - 1).toInt + entries(entries.length - 1).toInt - 1
  }
}

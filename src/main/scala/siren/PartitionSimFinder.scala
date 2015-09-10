package siren

import org.apache.spark.SparkContext
import SparkContext._

object PartitionSimFinder {
  def main(args : Array[String]): Unit = {
    args match {
      case Array("explorePartitionAdjacencies", clusterFile, genomeFile, readLenStr, minIntersectSizeStr, outFile) => {
        val readLen = readLenStr.toInt
        val minIntersectSize = minIntersectSizeStr.toInt

        // load genome
        println("Loading genome...")
        val genome = FASTA.read(genomeFile)
        val genomeSize = genome.totalSize
        println("Done loading genome.")
      
        // load cluster file
        // want pos => cluster id
        println("Loading cluster file...")
        val clusterIdByPos = ClusterHelper.getClusterIdByPos(clusterFile, true)
        println("Done loading cluster file.")
      
        // create partitions
        // create bitmap for genome
        // set any pos to true that is within a similar region
        println("Creating genome bitmap...")
        val (bitmapL, bitmapH) = getGenomeBitmaps(clusterFile, genomeSize, readLen)
        println("Done creating genome bitmap.")
      
        // 1st pass:  find out how many partitions there are
        println("Determining # of partitions...")
        val numPartitions = getNumberOfPartitions(bitmapL, bitmapH)
        println("# partitions: " + numPartitions)

        // 2nd pass:  for each partition, figure out which cluster IDs it has
        // represent each partition as a set of cluster IDs (use # partitions to allocate array of sets)
        println("Determining which cluster IDs make up each partition...")
        val partitions = getPartitionsByClusterIds(bitmapL, bitmapH, numPartitions, clusterIdByPos)
        println("Done with partitions and cluster IDs.")
      
        // do union find on the sets
        val unionFind = unionFindOnPartitions(partitions, minIntersectSize)
        
        // get clusters
        def isValid(arg1: Long, arg2: Long) = true
        unionFind.findClusters(0, isValid)
        
        // print clusters to file
        // want one line per cluster, largest to smallest
        // each line should have cluster size & then all partitions ids that belong to that cluster (should be sorted)
        val fw = new java.io.FileWriter(outFile)
        val bw = new java.io.BufferedWriter(fw)
        unionFind.writeClustersToFileSnapCompliant(bw, readLen, minIntersectSize, 2)
        bw.close
      }
      case Array("getPartitionClustersByPos", clusterFile, genomeFile, readLenStr, partitionClusterFile, outFile) => {
        val readLen = readLenStr.toInt

        // load genome
        println("Loading genome...")
        val genome = FASTA.read(genomeFile)
        val genomeSize = genome.totalSize
        println("Done loading genome.")
      
        // load cluster file
        // want pos => cluster id
        println("Loading cluster file...")
        val clusterIdByPos = ClusterHelper.getClusterIdByPos(clusterFile, true)
        println("Done loading cluster file.")
      
        // create partitions
        // create bitmap for genome
        // set any pos to true that is within a similar region
        println("Creating genome bitmap...")
        val (bitmapL, bitmapH) = getGenomeBitmaps(clusterFile, genomeSize, readLen)
        println("Done creating genome bitmap.")
      
        // 1st pass:  find out how many partitions there are
        println("Determining # of partitions...")
        val numPartitions = getNumberOfPartitions(bitmapL, bitmapH)
        println("# partitions: " + numPartitions)

        // 2nd pass:  for each partition, figure out which cluster IDs it has
        // represent each partition as a set of cluster IDs (use # partitions to allocate array of sets)
        println("Determining which cluster IDs make up each partition...")
        val partitions = getPartitionsByClusterIds(bitmapL, bitmapH, numPartitions, clusterIdByPos)
        println("Done with partitions and cluster IDs.")
        
        // For each partition cluster, give start/end positions for each partition belonging to that cluster
        // input format:  size partitionId1 partitionId2 ...
        // output format:  size (start1,end1) (start2,end2) ...
        // need partition ID => (startPos,endPos)
        val partitionEndPoints = getPartitionEndPoints(bitmapL, bitmapH, numPartitions)
        
        // Load partition clusters & open output file
        val source = scala.io.Source.fromFile(partitionClusterFile)
        val lines = source.getLines
        val fw = new java.io.FileWriter(outFile)
        val bw = new java.io.BufferedWriter(fw)
        
        var i = 0
        lines.foreach(l => {
          if (i == 0) {
            // print out header unchanged
            bw.write(l)
            bw.newLine
          } else {
            if ((i % 1000) == 0) println("Parsing line " + i + "...")
            val entries = l.split(" ")
            // first entry is size; subsequent entries are partition ids
            val size = entries(0)
            bw.write(size + " ")
            (1 until entries.length).foreach(n => {
              val partitionId = entries(n).toInt
              val endpoints = partitionEndPoints(partitionId)
              bw.write(endpoints.toString + " ")
            })
            bw.newLine
          }
          
          i += 1
        })
        
        bw.close
      }
      case Array("explorePartitionAdjacenciesParallel", sparkMaster, clusterFile, genomeFile, readLenStr, minIntersectSizeStr, outFile, gridDimStr) => {
        val readLen = readLenStr.toInt
        val minIntersectSize = minIntersectSizeStr.toInt
        val gridDim = gridDimStr.toInt

        // load genome
        println("Loading genome...")
        val genome = FASTA.read(genomeFile)
        val genomeSize = genome.totalSize
        println("Done loading genome.")
      
        // load cluster file
        // want pos => cluster id
        println("Loading cluster file...")
        val clusterIdByPos = ClusterHelper.getClusterIdByPos(clusterFile, true)
        println("Done loading cluster file.")
      
        // create partitions
        // create bitmap for genome
        // set any pos to true that is within a similar region
        println("Creating genome bitmap...")
        val (bitmapL, bitmapH) = getGenomeBitmaps(clusterFile, genomeSize, readLen)
        println("Done creating genome bitmap.")
      
        // 1st pass:  find out how many partitions there are
        println("Determining # of partitions...")
        val numPartitions = getNumberOfPartitions(bitmapL, bitmapH)
        println("# partitions: " + numPartitions)

        // 2nd pass:  for each partition, figure out which cluster IDs it has
        // represent each partition as a set of cluster IDs (use # partitions to allocate array of sets)
        println("Determining which cluster IDs make up each partition...")
        val partitions = getPartitionsByClusterIds(bitmapL, bitmapH, numPartitions, clusterIdByPos)
        println("Done with partitions and cluster IDs.")
      
        // Do union find
        // BEGIN SPARK PART
        println("Begin spark part...")
        
        // Divide up partitions into grid cells
        //val numPartitions = 650000 // DEBUG
        val rangeLength = (numPartitions / gridDim).toInt
        val rangeStarts = (0 until numPartitions by rangeLength)
        println("Ranges:")
        rangeStarts.foreach(rangeStart => {
          val startPartition = rangeStart
          val endPartition = math.min(rangeStart + rangeLength, numPartitions) - 1
          println(startPartition + ", " + endPartition)
        })
        val gridStarts = rangeStarts.map(i => List(i).padTo(rangeStarts.length, i).zip(rangeStarts)).flatten  // gives start for indexing & scanning range

        // Initialize union find (using UnionFindL)
        var ufClusters = new UnionFindL(numPartitions) // must use "L" here b/c it works with accumulator classes
        
        // Initialize spark context & accumulator
        val sc = new SparkContext(sparkMaster, "PartitionSimFinder", "", Seq("target/scala-2.9.2/snap_2.9.2-0.0.jar"))
        val ufAccumulator = sc.accumulator(ufClusters.asInstanceOf[UnionFindAbstract])(UnionFindAP)
        
        // Parallelize -- requires function to do union find on a grid cell
        sc.parallelize(gridStarts, gridStarts.size).foreach(p => {
          println("In task " + p + "...")
          val uf = getGridCellClusters(p, rangeLength, partitions, minIntersectSize)
          ufAccumulator += uf
        })
        
        // END SPARK PART
        println("End spark part...")
        
        // get clusters
        println("Finding clusters...")
        ufClusters = ufAccumulator.value.asInstanceOf[UnionFindL]
        def isInvalid(arg: Array[Byte]) = false // hack
        ufClusters.findClusters(0, isInvalid) // not sure 0 will work for readLen here
        println("Done finding clusters.")
        
        // print clusters to file
        // want one line per cluster, largest to smallest
        // each line should have cluster size & then all partitions ids that belong to that cluster (should be sorted)
        println("Printing clusters to file...")
        val fw = new java.io.FileWriter(outFile)
        val bw = new java.io.BufferedWriter(fw)
        // unionFind.writeClustersToFileSnapCompliant(bw, readLen, minIntersectSize, 2)
        // From SimFinder: ufClusters.writeClustersToFileSnapCompliant(bw, GenomeLoader.genome, params.readLen, params.unionDist, params.minClusterSize)
        ufClusters.writeClustersToFileSnapCompliant(bw, genome, readLen, minIntersectSize, 2 /* minClusterSize */)
        bw.close
        println("Done printing.")
      }
      case Array("explorePartitionAdjacenciesParallelV2", sparkMaster, clusterFile, genomeFile, readLenStr, minIntersectSizeStr, outFile, gridDimStr) => {
        val readLen = readLenStr.toInt
        val minIntersectSize = minIntersectSizeStr.toInt
        val gridDim = gridDimStr.toInt

        // load genome
        println("Loading genome...")
        val genome = FASTA.read(genomeFile)
        val genomeSize = genome.totalSize
        println("Done loading genome.")
        
      }
      case Array("testSparkInit") => {
        val sc = new SparkContext("local[16]", "PartitionSimFinderTester", "", Seq("target/scala-2.9.2/snap_2.9.2-0.0.jar"))
        
        // Parallelize -- requires function to do union find on a grid cell
        sc.parallelize(1 to 10).foreach(p => {
          println("In task " + p + "...")
        })
      }
      case Array("formatConverter", genomeFile, inFile, outFile) => {
        // Load genome
        val genome = FASTA.read(genomeFile)

        // Open input & output files
        val source = scala.io.Source.fromFile(inFile)
        val lines = source.getLines
        val fw = new java.io.FileWriter(outFile)
        val bw = new java.io.BufferedWriter(fw)
        
        // Parse input file & write newly-formatted data to output
        val endPoints = """\((\d+),(\d+)\)""".r
        var i = 0
        lines.foreach(l => {
          // skip header
          if (i > 0) {
            val entries = l.split(" ")
            // skip size
            val familySize = entries(0).toInt
            (1 until entries.length).foreach(j => {
              val endPoints(start, end) = entries(j)

              val (startPiece, startPos) = genome.getLocation(start.toLong)
              val (endPiece, endPos) = genome.getLocation(end.toLong)
              
              assert(startPiece == endPiece)
              
              bw.write("(" + List(startPiece, startPos, endPos).mkString(", ") + ")" + " ")
            })
            bw.newLine
          }
          i += 1
        })

        bw.close
      }
      case _ => println("Incorrect args.")
    }
  }
  
  def getNumberOfPartitions(bitmapL: Array[Int], bitmapH: Array[Int]): Int = {
    var numPartitions = 0
    var prevPos = 0
    var i = 0L
    val midPoint = getMidpointFromBitmaps(bitmapL, bitmapH)
    val genomeSize = getGenomeSizeFromBitmaps(bitmapL, bitmapH)
    
    while (i < genomeSize) {
      val bitmapVal = 
       if (i < midPoint) bitmapL(i.toInt)
       else bitmapH((i - midPoint).toInt)

      // if prevPos == bitmapVal == 0, do nothing
      if (prevPos == 0 && bitmapVal == 1) {
        // you've reached the beginning of a similar region
        //numPartitions += 1 // NOTE:  I think this is wrong... should be in next if clause; may also need to be updated in analyzeContiguousSimilarRegions
        prevPos = 1
      }
      else if (prevPos == 1 && (bitmapVal == 0 || i == (genomeSize - 1))) {
        // you've reached the end of a similar region OR the end of the genome
        prevPos = 0
        numPartitions += 1
      }
      // if prevPos == bitmapVal == 1, do nothing

      i += 1
    }

    numPartitions
  }

  def getPartitionsByClusterIds(bitmapL: Array[Int], bitmapH: Array[Int], numPartitions: Int, clusterIdByPos: scala.collection.mutable.Map[Long, Long]): Array[scala.collection.mutable.Set[Long]] = {
    val partitions = Array.fill(numPartitions)(scala.collection.mutable.Set[Long]())
    var partitionStartPos = 0L
    var partitionNum = 0
    
    val midPoint = getMidpointFromBitmaps(bitmapL, bitmapH)
    val genomeSize = getGenomeSizeFromBitmaps(bitmapL, bitmapH)
    
    var i = 0L
    var prevPos = 0
    while (i < genomeSize) {
      val bitmapVal = 
        if (i < midPoint) bitmapL(i.toInt)
        else bitmapH((i - midPoint).toInt)

      // if prevPos == bitmapVal == 0, do nothing
      if (prevPos == 0 && bitmapVal == 1) {
        // you've reached the beginning of a similar region OR the end of the genome
        prevPos = 1
        partitionStartPos = i
      }
      else if (prevPos == 1 && (bitmapVal == 0 || i == (genomeSize - 1))) {
        // you've reached the end of a similar region OR the end of the genome
        prevPos = 0
        
        // we have come to the end of a partition; now, add it to the partitions array
        //partitions(partitionNum) = scala.collection.mutable.Set[Long]()
        var k = partitionStartPos
        while (k < i) {
          // check whether this pos is in a similar region
          // if so, add its id to this partition's set
          val id = clusterIdByPos.get(k)
          if (id != None) {
            partitions(partitionNum) += id.get
          }
          
          k += 1
        }
        
        partitionNum += 1
      }
      // if prevPos == bitmapVal == 1, do nothing

      i += 1
    }
    
    partitions
  }

  def getGenomeSizeFromBitmaps(bitmapL: Array[Int], bitmapH: Array[Int]): Long = bitmapL.length.toLong + bitmapH.length.toLong
  
  def getMidpointFromBitmaps(bitmapL: Array[Int], bitmapH: Array[Int]): Int = bitmapL.length

  def getMidPointFromGenomeSize(genomeSize: Long): Int = (genomeSize / 2).toInt
  
  def getGenomeBitmaps(clusterFile: String, genomeSize: Long, readLen: Int): (Array[Int], Array[Int]) = {
    val midPoint = getMidPointFromGenomeSize(genomeSize)

    // 0 => not in similar region
    // 1 => in similar region
    // keep everything 0-indexed (both genome and similar regions are)
    val bitmapL = Array.fill(midPoint)(0) // 0 to midPoint - 1
    val bitmapH = Array.fill((genomeSize - midPoint).toInt)(0) // midPoint to genomeSize - 1
    
    // read in cluster file
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
    
    (bitmapL, bitmapH)
  }

  def getPartitionIdsByClusterIds(partitions: Array[scala.collection.mutable.Set[Long]]): scala.collection.mutable.Map[Long, List[Int]] = {
    getPartitionIdsByClusterIds(partitions, 0, partitions.length)
  }
  
  def getPartitionIdsByClusterIds(partitions: Array[scala.collection.mutable.Set[Long]], start: Int /* inclusive */, end: Int /* exclusive */): scala.collection.mutable.Map[Long, List[Int]] = {
    assert(start >= 0)
    assert(end <= partitions.length)
    val partitionIdsByClusterId = scala.collection.mutable.Map[Long, List[Int]]()
    (start until end).foreach(p => {
      val ids = partitions(p)
      ids.foreach(i => {
        partitionIdsByClusterId.get(i) match {
          case None => {
            // haven't yet encountered this cluster id => make list for it
            partitionIdsByClusterId += ((i, List(p)))
          }
          case Some(list) => {
            // this cluster id already has a list => get it, and append to it
            val newList = p :: list
            partitionIdsByClusterId += ((i, newList))
          }
        }
      })
    })
    
    partitionIdsByClusterId
  }

  def unionFindOnPartitions(partitions: Array[scala.collection.mutable.Set[Long]], minIntersectSize: Int): UnionFind = {
    // use partitionIdsByClusterId to avoid doing numPartitions x numPartitions similarity ops
    val partitionIdsByClusterId = getPartitionIdsByClusterIds(partitions)
    val numPartitions = partitions.length
    val unionFind = new UnionFind(numPartitions)
    
    (0 until numPartitions).foreach(p => {
      if (p % 10000 == 0)
        println("On partition " + p + "...")
      // for each cluster id in this partition, get the other partitions that have this cluster id
      // compare current partition to each partition from map
      // if intersection is sufficiently large, merge
      val clusterIds = partitions(p)
      clusterIds.foreach(i => {
        val otherPartitionIds = partitionIdsByClusterId(i)  // this cluster id should be there, so if not, crash is okay
        otherPartitionIds.foreach(otherP => {
          if (clusterIds.intersect(partitions(otherP)).size >= minIntersectSize)
            unionFind.union(p, otherP)
        })
      })
    })

    unionFind
  }

  def getPartitionEndPoints(bitmapL: Array[Int], bitmapH: Array[Int], numPartitions: Int): Array[(Long, Long)] = {
    val partitionEndPoints = Array.fill(numPartitions)((0L, 0L))
    var partitionNum = 0
    var prevPos = 0
    var partitionStartPos = 0L
    var i = 0L
    val midPoint = getMidpointFromBitmaps(bitmapL, bitmapH)
    val genomeSize = getGenomeSizeFromBitmaps(bitmapL, bitmapH)
    
    while (i < genomeSize) {
     val bitmapVal = 
       if (i < midPoint) bitmapL(i.toInt)
       else bitmapH((i - midPoint).toInt)

     // if prevPos == bitmapVal == 0, do nothing
     if (prevPos == 0 && bitmapVal == 1) {
       // you've reached the beginning of a partition
       partitionStartPos = i
       prevPos = 1
     }
     else if (prevPos == 1 && (bitmapVal == 0 || i == (genomeSize - 1))) {
       // you've reached the end of a partition OR the end of the genome
       // => store this partition
       var endPos = 
         if (i == (genomeSize - 1)) i
         else i - 1
       partitionEndPoints(partitionNum) = ((partitionStartPos, endPos))
       
       partitionNum += 1
       prevPos = 0
     }
     // if prevPos == bitmapVal == 1, do nothing

     i += 1
    }

    partitionEndPoints
  }

  def getGridCellClusters(gridStart: (Int, Int), rangeLength: Int, partitions: Array[scala.collection.mutable.Set[Long]], minIntersectSize: Int): UnionFindAbstract = {
    val numPartitions = partitions.length

    // Figure out which grid cell I'm on
    val indexStart = gridStart._1  // range start for row, ie, which positions to index // inclusive
    val indexEnd = math.min(indexStart + rangeLength, numPartitions)  // exclusive

    val scanStart = gridStart._2   // range start for column, ie, which positions to scan over for clustering // inclusive
    val scanEnd = math.min(scanStart + rangeLength, numPartitions)  // exclusive

    println("Indexing " + indexStart + " to " + (indexEnd - 1) + "; scanning " + scanStart + " to " + (scanEnd - 1) + "...")

    // Build index for given range
    val partitionIdsByClusterId = getPartitionIdsByClusterIds(partitions, indexStart, indexEnd)

    // Init
    val uf = 
      if (indexStart == scanStart) new UnionFindGridDiagonal((indexStart, indexEnd - 1 /* expects inclusive */))
      else new UnionFindGrid((indexStart, indexEnd), (scanStart, scanEnd))
    
    (scanStart until scanEnd).foreach(p => {
      // for each cluster id in this partition, get the other partitions that have this cluster id
      // compare current partition to each partition from map
      // if intersection is sufficiently large, merge
      val clusterIds = partitions(p)
      clusterIds.foreach(i => {
        val otherPartitionIds = partitionIdsByClusterId.get(i)
        otherPartitionIds match {
          case Some(ids) => 
            ids.foreach(otherP => {
              if (clusterIds.intersect(partitions(otherP)).size >= minIntersectSize)
                uf.union(p, otherP)
            })
          case None =>
        }
      })
    })

    def isInvalid(arg: Array[Byte]) = false // hack
    uf.findClusters(0, isInvalid) // not sure 0 will work for readLen here
    uf
  }
}

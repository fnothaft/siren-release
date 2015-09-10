package siren

import it.unimi.dsi.fastutil.longs.LongList
import it.unimi.dsi.fastutil.longs.LongArrayList

class UnionFind(numBases: Long) {
  val parent = (0 until numBases.toInt).toArray
  val rank = Array.fill(numBases.toInt)(0)
  val clusterSize = Array.fill(numBases.toInt)(0) // how big is the cluster that this pos belongs to?
  var totalClusters = 0
  var nonTrivialClusters = List[Int]()
  var nonTrivialMembers = scala.collection.mutable.Map[Long, LongArrayList]()
  var firstMember = Array.fill(numBases.toInt)(0)

  def find(v: Int): Int = {
    if (v != parent(v)) {
      parent(v) = find(parent(v))
    }
    return parent(v)
  }
  
  def union(v: Int, w: Int) = {
    val x = find(v)
    val y = find(w)

    if (rank(x) > rank(y)) {
      parent(y) = x
    } else {
      parent(x) = y
      if (rank(y) == rank(x))
        rank(y) += 1
    }
  }
  
  // assumes nontrivial cluster size is 2
  def findClusters(readLen: Int, isValid: ((Long, Long) => Boolean)) = {
    var pos = 0
    while (pos < numBases - readLen) {
      if (isValid(pos, pos + readLen)) {
        val p = find(pos)
        clusterSize(p) += 1
        
        if (clusterSize(p) == 1) {
          totalClusters += 1
          firstMember(p) = pos  // Remember in case cluster becomes non-trivial
        } else if (clusterSize(p) == 2) {
          // Just found a new non-trivial cluster (with more than one element)
          nonTrivialClusters = p :: nonTrivialClusters
          val l = new LongArrayList
          l.add(firstMember(p))
          l.add(pos)
          nonTrivialMembers += ((p, l))
        } else if (clusterSize(p) > 2) {
          val members = nonTrivialMembers.get(p).get
          members.add(pos)
        }
      }
      
      pos += 1
    }
    
    nonTrivialClusters = nonTrivialClusters.sortWith((c1, c2) => clusterSize(c1) > clusterSize(c2))
  }

  def printClusterStats(validPositions: Int) = {
    println("Valid positions: " + validPositions)
    println("Total clusters: " + totalClusters)
    println("Non-trivial clusters: " + nonTrivialClusters.size)

    val nonTrivialPositions = validPositions - (totalClusters - nonTrivialClusters.size)
    
    println("Positions in non-trivial clusters: " + nonTrivialPositions)
    println("Mean size of non-trivial clusters: " + nonTrivialPositions / nonTrivialClusters.size.toDouble)
  }
  
  def writeClustersToFile(bw: java.io.BufferedWriter, genome: Genome, readLen: Int, minClusterSize: Int = 5) = {
    var numPrinted = 0
    var i = 0
    val clusterRep = new Array[Byte](readLen)
    val buf1 = new Array[Byte](readLen)
    val buf2 = new Array[Byte](readLen)
    
    while (i < nonTrivialClusters.size) {
      val id = nonTrivialClusters(i)
      if (clusterSize(id) >= minClusterSize) {
        genome.getSubstring(id, id + readLen, clusterRep)  // get first member of cluster; should eventually be consensus string
        val members = nonTrivialMembers(id)
        // Compute maximum & average Hamming distance among pairs in cluster
        var maxDist = 0
        var sumDists = 0
        var j = 0
        while (j < members.size) {
          genome.getSubstring(members.getLong(j), members.getLong(j) + readLen, buf1)
          var k = j + 1
          while (k < members.size) {
            genome.getSubstring(members.getLong(k), members.getLong(k) + readLen, buf2)
            val dist = HammingDistance.distance(buf1, buf2, readLen)
            maxDist = math.max(dist, maxDist)
            sumDists += dist
            
            k += 1
          }
          
          j += 1
        }
        
        val avgDist = sumDists.toDouble / (members.size * (members.size - 1) / 2)

        // Print info about the cluster
        bw.write(List(members.size, avgDist, maxDist, readLen, new String(clusterRep)).mkString(" "))
        j = 0
        while (j < members.size) {
          bw.write(" " + members.getLong(j))
          
          j += 1
        }
        bw.newLine
        
        numPrinted += 1
      } else {
        i = nonTrivialClusters.size
      }
      
      i += 1
    }

    println("Printed out " + numPrinted + " clusters of size > " + minClusterSize)
  }

  def writeKernelMatrixToFile(bw: java.io.BufferedWriter, genome: Genome, readLen: Int, members: LongArrayList) = {
    val buf1 = new Array[Byte](readLen)
    val buf2 = new Array[Byte](readLen)
    
    // Compute Hamming distance for each pair in cluster
    var j = 0
    while (j < members.size) {
      genome.getSubstring(members.getLong(j), members.getLong(j) + readLen, buf1)
      var k = 0
      while (k < members.size) {
        genome.getSubstring(members.getLong(k), members.getLong(k) + readLen, buf2)
        val dist = HammingDistance.distance(buf1, buf2, readLen)
        bw.write(" " + dist)
        
        k += 1
      }
      bw.newLine
      
      j += 1
    }
  }
  
  def writeTopNKernelMatricesToFile(fnamePrefix: String, genome: Genome, readLen: Int, n: Int) = {
    (0 until n).foreach(c => {
      val outKernel = new java.io.File(fnamePrefix + ".kernel." + c)
      if (!outKernel.exists)
        outKernel.createNewFile

      val fwKernel = new java.io.FileWriter(outKernel.getName)
      val bwKernel = new java.io.BufferedWriter(fwKernel)
      
      val clusterId = nonTrivialClusters(c)
      val members = nonTrivialMembers(clusterId)
      writeKernelMatrixToFile(bwKernel, genome, readLen, members)
      
      bwKernel.close
    })
  }

  def writeClusterMembersToFile(bw: java.io.BufferedWriter, genome: Genome, readLen: Int, members: LongArrayList) = {
    val buf = new Array[Byte](readLen)
    
    (0 until members.size).foreach(i => {
      genome.getSubstring(members.getLong(i), members.getLong(i) + readLen, buf)
      bw.write(i + " " + (new String(buf)))
      bw.newLine
    })
  }

  def writeTopNClusterMembersToFile(fnamePrefix: String, genome: Genome, readLen: Int, n: Int) = {
    (0 until n).foreach(c => {
      val out = new java.io.File(fnamePrefix + ".members." + c)
      if (!out.exists)
        out.createNewFile

      val fw = new java.io.FileWriter(out.getName)
      val bw = new java.io.BufferedWriter(fw)
      
      val clusterId = nonTrivialClusters(c)
      val members = nonTrivialMembers(clusterId)
      writeClusterMembersToFile(bw, genome, readLen, members)
      
      bw.close
    })
  }

  // output format expected by SNAP
  // Line 1:  <genomeLength> <numClusters> <stringLength> <mergeDist>
  // For each cluster:  <numMembers> <member1> <member2> <...>
  def writeClustersToFileSnapCompliant(bw: java.io.BufferedWriter, readLen: Int, mergeDist: Int, minClusterSize: Int = 5 /* inclusive */) = {
    var numPrinted = 0
    var i = 0
    
    // print header line
    // numClusters should be # clusters that I'm printing out (ie, with size above minClusterSize)
    bw.write(List(numBases, nonTrivialClusters.size, readLen, mergeDist).mkString(" "))
    bw.newLine

    while (i < nonTrivialClusters.size) {
      val id = nonTrivialClusters(i)
      val clusterSizeId = clusterSize(id)
      
      if (clusterSizeId >= minClusterSize) {
        val members = nonTrivialMembers(id)

        // Print info about the cluster
        //bw.write(List(members.size, avgDist, maxDist, readLen, new String(clusterRep)).mkString(" "))
        bw.write(members.size.toString)
        var j = 0
        while (j < members.size) {
          bw.write(" " + members.getLong(j))
          j += 1
        }
        bw.newLine
        numPrinted += 1
      } else {
        i = nonTrivialClusters.size
      }
      
      i += 1
    }

    println("Printed out " + numPrinted + " clusters of size >= " + minClusterSize)
  }
}

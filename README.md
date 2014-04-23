#Modularity Clustering#
A version of modularity clustering algorithm by Newman Phys Rev E 69, 066133 PDF Here: 
http://arxiv.org/pdf/cond-mat/0309508.pdf

To run it via python import it as a module (you'll want to make sure the correct folder is in your ```sys.path```):
```python
from modularityClustering import modCluster
```

See file:
```bash
examples/exampleModularityClustering.py
```
For more information see comments in ```examples/exmapleModularityClustering.py``` .

Output will be of three types: 
  1.  A JSON file which includes the different clusters and the modularity values for each cluster. 
  ```json 
  {
      "0": {
          "a": 0.3, 
          "e": {
              "0": 0.1, 
              "1": 0.0001
          }, 
          "members": [
              "myNodeName1", 
              "myNodeName2"
          ]
      }, 
      "1": {
          "a": 0.4, 
          "e": {
              "0": 0.001, 
              "1": 0.3
          }, 
          "members": [
              "myNodeName3", 
              "myNodeName4",
              "myNodeName5"
          ]
      }
  }
  ```

  2.  A five column filename where you can make a pretty gnuplot correlation plot.
      To plot, open gnuplot and use following commands:
      
      ```bash
      $set pm3d map; unset key
      $splot "myDataOutput.dat" u 3:4:5
      $set terminal pngcairo size 640,640 enhanced font 'Verdana,10'
      $set border 3 back; set bmargin 10; set lmargin 10; set border
      $set output "myOutputPlot.png"
      $replot; unset out; exit
      ```

  3.  And a file useful for making a scala map which takes the feature and points to it's cluster ID with the following format ```bash "stringID" -> clusterID``` -- ie:
  
  ```bash
  "myNodeName" -> 0
  "otherNodeName" -> 0
  "yetAnotherNodeName" -> 1
  ```

You can run it on the data provided (3 Column CSV of format Node,Node,EdgeWeight)
```bash
karate.csv
```
This is the famous "Zachary's Karate Club" data from 1977.

Here's a good example of how to use the module:

```python
from modularityClustering import modCluster
mc = modCluster()
mc.loadEdges("karate.csv")
mc.findCommunities()
mc.printClustersJSON("karate.json")
```


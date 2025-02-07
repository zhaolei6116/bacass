nextflow.enable.dsl=2

def getAbsPath(String filename) {
  return new File("$PWD", filename).absolutePath
}
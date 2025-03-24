nextflow.enable.dsl=2
import groovy.json.JsonSlurper


def getAbsPath(String filename) {
  return new File("$PWD", filename).absolutePath
}


// D定义全局变量



def getBarcode(barcode_version, barcode_id) {
  def jsonSlurper = new JsonSlurper()
  def barcodes = jsonSlurper.parse(new File('/nas02/project/zhaolei/pipeline/bacteria_genome_assembly/lib/NBD114-96.json'))
  // check
  if (!barcode_version){
    throw new Exception("Barcode version ${barcode_version} not found in barcodes.json")
  }
  // check
  if (!barcode_id){
    throw new Exception("Barcode ID ${barcode_id} not found in barcodes.json")
  }
  // get back seq
  return barcodes[barcode_version][barcode_id]
}
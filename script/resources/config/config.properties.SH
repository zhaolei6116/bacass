url=https://www.cwbio.com/gwapi/get_report
responseUrl=https://www.cwbio.com/gwapi/push_report
appid=cwbioprimecx
sign=5dbc98dcc983a70728bd082d1a47546e
appsecret=5dbc98dcc983a70728bd082d1a47546e

startTime=2024-12-1 08:00:00
endTime=2024-12-29 20:00:00
bufferSize=8192
batchSize=100
rateLimit=100.0

maxRetries=5
initialDelayMs=1000
backoffMultiplier=2.0
maxDelayMs=30000
retryDelaySeconds=10
timeoutSeconds=30

downloadPath=downloads
uploadPath=uploads
tempPath=temp

metrics.enabled=true
metrics.reportIntervalSeconds=60

logging.level.root=INFO

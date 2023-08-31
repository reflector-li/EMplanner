#include "server_config.h"

namespace foxglove {
template <>
void Server<WebSocketNoTls>::setupTlsHandler() {}
}  // namespace foxglove

std::string FoxgloveServer::SerializeFdSet(const google::protobuf::Descriptor* toplevelDescriptor) {
  google::protobuf::FileDescriptorSet fdSet;
  std::queue<const google::protobuf::FileDescriptor*> toAdd;
  toAdd.push(toplevelDescriptor->file());
  std::unordered_set<std::string> seenDependencies;
  while (!toAdd.empty()) {
    const google::protobuf::FileDescriptor* next = toAdd.front();
    toAdd.pop();
    next->CopyTo(fdSet.add_file());
    for (int i = 0; i < next->dependency_count(); ++i) {
      const auto& dep = next->dependency(i);
      if (seenDependencies.find(dep->name()) == seenDependencies.end()) {
        seenDependencies.insert(dep->name());
        toAdd.push(dep);
      }
    }
  }
  return fdSet.SerializeAsString();
}



FoxgloveServer::FoxgloveServer(){
  const auto logHandler = [](foxglove::WebSocketLogLevel, char const* msg) {
    std::cout << msg << std::endl;
  };
  foxglove::ServerOptions serverOptions;
  server_ = std::make_unique<foxglove::Server<foxglove::WebSocketNoTls>>(
    "C++ Protobuf example server", logHandler, serverOptions);

  foxglove::ServerHandlers<foxglove::ConnHandle> hdlrs;
  hdlrs.subscribeHandler = [&](foxglove::ChannelId chanId, foxglove::ConnHandle) {
    std::cout << "first client subscribed to " << chanId << std::endl;
  };
  hdlrs.unsubscribeHandler = [&](foxglove::ChannelId chanId, foxglove::ConnHandle) {
    std::cout << "last client unsubscribed from " << chanId << std::endl;
  };
  server_->setHandlers(std::move(hdlrs));
}


void FoxgloveServer::start(const std::string& host, uint16_t port){
    server_->start(host, port);
}


void FoxgloveServer::stop(){
    server_->stop();
}

std::vector<foxglove::ChannelId> FoxgloveServer::addChannels(const std::vector<foxglove::ChannelWithoutId>& channels){
   std::vector<foxglove::ChannelId> ChannelIds = server_->addChannels({channels});
   return ChannelIds;
}

void FoxgloveServer::removeChannels(const std::vector<foxglove::ChannelId>& ChannelIds){
    server_->removeChannels({ChannelIds});
}

void FoxgloveServer::sendMessages(foxglove::ChannelId chanId, uint64_t timestamp,const uint8_t* payload,size_t payloadSize){
  server_->broadcastMessage(chanId, timestamp, payload, payloadSize);
}

std::unique_ptr<foxglove::Server<foxglove::WebSocketNoTls>> & FoxgloveServer::getWebServer(){
  return server_;
}
#ifndef SERVER_CONFIG_H
#define SERVER_CONFIG_H

#include <string>
#include <vector>
#include <foxglove/websocket/base64.hpp>
#include <foxglove/websocket/common.hpp>
#include <foxglove/websocket/websocket_notls.hpp>
#include <foxglove/websocket/websocket_server.hpp>
#include <memory>

#include <google/protobuf/descriptor.pb.h>
#include <google/protobuf/util/time_util.h>

class FoxgloveServer
{
public:
    std::unique_ptr<foxglove::Server<foxglove::WebSocketNoTls>> server_;


    FoxgloveServer();
    ~FoxgloveServer() = default;
    void start(const std::string& host, uint16_t port);
    void stop();
    std::vector<foxglove::ChannelId> addChannels(const std::vector<foxglove::ChannelWithoutId>& channels);
    void removeChannels(const std::vector<foxglove::ChannelId>& ChannelIds);
    void sendMessages(foxglove::ChannelId chanId, uint64_t timestamp,const uint8_t* payload,size_t payloadSize);
    std::unique_ptr<foxglove::Server<foxglove::WebSocketNoTls>> & getWebServer();

    static std::string SerializeFdSet(const google::protobuf::Descriptor* toplevelDescriptor);
    
};

#endif

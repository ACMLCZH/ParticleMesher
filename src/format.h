#pragma once

#include <string>
#include <vector>
#include <numeric>

namespace particle_mesher {

template<typename ... Args>
static std::string str_format(const std::string& format, Args ... args) {
	auto size_buf = std::snprintf(nullptr, 0, format.c_str(), args ...) + 1;
	std::unique_ptr<char[]> buf(new(std::nothrow) char[size_buf]);

	if (!buf)
		return std::string("");

	std::snprintf(buf.get(), size_buf, format.c_str(), args ...);
	return std::string(buf.get(), buf.get() + size_buf - 1);
}

static void combine_message(const std::vector<std::string> &msg_pack, std::string &msg_str) noexcept {
    auto tot_len = std::accumulate(
        msg_pack.begin(), msg_pack.end(), 0u,
        [](const uint &prev, const std::string &info) { return prev + info.size(); } 
    );
    msg_str.clear();
    msg_str.reserve(tot_len + 1);
    for (const auto &msg: msg_pack) msg_str += msg;
}

}// namespace particle_mesher


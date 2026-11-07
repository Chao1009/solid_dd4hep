// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// Minimal SHA-256 (FIPS 180-4), used to check downloaded data files (field maps) against the checksum given in the
// compact file. No external dependency.
//   solid::SHA256 h; h.update(data, n); std::string hex = h.hexdigest();   // 64 lower-case hex digits
//   std::string hex = solid::FileSHA256("fieldmaps/map.dat");               // throws if the file cannot be read

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace solid {

class SHA256 {
public:
  SHA256() { reset(); }

  void reset()
  {
    m_h    = {0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au, 0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u};
    m_len  = 0;
    m_nbuf = 0;
  }

  void update(const unsigned char* data, std::size_t n)
  {
    m_len += n;
    while (n > 0) {
      const std::size_t take = std::min<std::size_t>(64 - m_nbuf, n);
      std::memcpy(m_buf.data() + m_nbuf, data, take);
      m_nbuf += take;
      data += take;
      n -= take;
      if (m_nbuf == 64) {
        block(m_buf.data());
        m_nbuf = 0;
      }
    }
  }

  // finishes the hash (the object must be reset() before it is used again)
  std::string hexdigest()
  {
    const std::uint64_t        bits = m_len * 8;
    static const unsigned char zero = 0x00;
    const unsigned char        one  = 0x80;
    update(&one, 1);
    while (m_nbuf != 56)
      update(&zero, 1);
    unsigned char len_be[8];
    for (int i = 0; i < 8; ++i)
      len_be[i] = static_cast<unsigned char>(bits >> (56 - 8 * i));
    update(len_be, 8);
    static const char hex[] = "0123456789abcdef";
    std::string       out;
    out.reserve(64);
    for (std::uint32_t v : m_h)
      for (int s = 28; s >= 0; s -= 4)
        out += hex[(v >> s) & 0xfu];
    return out;
  }

private:
  static std::uint32_t rotr(std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

  void block(const unsigned char* p)
  {
    static const std::uint32_t k[64] = {
        0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
        0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
        0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
        0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
        0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
        0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
        0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
        0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};
    std::uint32_t w[64];
    for (int i = 0; i < 16; ++i)
      w[i] = (std::uint32_t(p[4 * i]) << 24) | (std::uint32_t(p[4 * i + 1]) << 16) |
             (std::uint32_t(p[4 * i + 2]) << 8) | std::uint32_t(p[4 * i + 3]);
    for (int i = 16; i < 64; ++i) {
      const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
      const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i]                   = w[i - 16] + s0 + w[i - 7] + s1;
    }
    std::uint32_t a = m_h[0], b = m_h[1], c = m_h[2], d = m_h[3], e = m_h[4], f = m_h[5], g = m_h[6], h = m_h[7];
    for (int i = 0; i < 64; ++i) {
      const std::uint32_t t1 = h + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + k[i] + w[i];
      const std::uint32_t t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
      h                      = g;
      g                      = f;
      f                      = e;
      e                      = d + t1;
      d                      = c;
      c                      = b;
      b                      = a;
      a                      = t1 + t2;
    }
    m_h[0] += a;
    m_h[1] += b;
    m_h[2] += c;
    m_h[3] += d;
    m_h[4] += e;
    m_h[5] += f;
    m_h[6] += g;
    m_h[7] += h;
  }

  std::array<std::uint32_t, 8>  m_h{};
  std::array<unsigned char, 64> m_buf{};
  std::uint64_t                 m_len  = 0;
  std::size_t                   m_nbuf = 0;
};

// SHA-256 of a file's contents as 64 lower-case hex digits
inline std::string FileSHA256(const std::string& path)
{
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw std::runtime_error("cannot read " + path);
  SHA256            h;
  std::vector<char> buf(1 << 20);
  while (in) {
    in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
    const std::streamsize n = in.gcount();
    if (n > 0)
      h.update(reinterpret_cast<const unsigned char*>(buf.data()), static_cast<std::size_t>(n));
  }
  if (in.bad())
    throw std::runtime_error("read error in " + path);
  return h.hexdigest();
}

} // namespace solid

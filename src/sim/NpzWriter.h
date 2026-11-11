// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 SoLID collaboration
//
// Minimal writer of numpy .npz files: a zip archive with one .npy member per array, readable with numpy.load().
// Members are deflated with zlib when SOLID_HAVE_ZLIB is defined and stored otherwise; both are valid .npz.
// Limits: little-endian hosts, archives and members below 4 GB (no zip64), at most 65535 members.
//
//   solid::NpzWriter npz;
//   npz.add("node_matrix", matrices /* std::vector<double> */, {n, 12});
//   npz.add("node_name", names /* std::vector<std::string>, UTF-8 */);
//   npz.addScalar("meta", json);
//   npz.write("scene.npz");          // throws std::runtime_error on failure

#pragma once

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#ifdef SOLID_HAVE_ZLIB
#include <zlib.h>
#endif

namespace solid {

class NpzWriter {
public:
  // numeric array; shape defaults to 1-D
  template <typename T> void add(const std::string& name, const std::vector<T>& data, std::vector<std::size_t> shape = {})
  {
    static_assert(std::is_arithmetic_v<T> && !std::is_same_v<T, bool>, "use std::vector<std::uint8_t> with addBool()");
    if (shape.empty()) {
      shape = {data.size()};
    }
    addRaw(name, descr<T>(), shape, data.data(), data.size() * sizeof(T));
  }

  // boolean array from bytes (0 or 1)
  void addBool(const std::string& name, const std::vector<std::uint8_t>& data)
  {
    addRaw(name, "|b1", {data.size()}, data.data(), data.size());
  }

  // 1-D array of strings (UTF-8 in, numpy fixed-width unicode '<U<n>' out)
  void add(const std::string& name, const std::vector<std::string>& data)
  {
    std::vector<std::u32string> u;
    std::size_t                 width = 1;
    for (const auto& s : data) {
      u.push_back(decodeUtf8(s));
      width = std::max(width, u.back().size());
    }
    std::vector<char32_t> buf(width * data.size(), U'\0');
    for (std::size_t i = 0; i < u.size(); ++i) {
      std::copy(u[i].begin(), u[i].end(), buf.begin() + i * width);
    }
    addRaw(name, "<U" + std::to_string(width), {data.size()}, buf.data(), buf.size() * sizeof(char32_t));
  }

  // 0-dimensional string (e.g. a JSON document)
  void addScalar(const std::string& name, const std::string& value)
  {
    std::u32string u = decodeUtf8(value);
    if (u.empty()) {
      u.push_back(U'\0');
    }
    addRaw(name, "<U" + std::to_string(u.size()), {}, u.data(), u.size() * sizeof(char32_t), true);
  }

  void write(const std::string& path) const
  {
    if (m_members.size() > 0xffff) {
      throw std::runtime_error("NpzWriter: too many members");
    }
    std::string out;
    std::string central;
    for (const auto& m : m_members) {
      if (out.size() + m.packed.size() > 0xfff00000u) {
        throw std::runtime_error("NpzWriter: archive larger than 4 GB");
      }
      const std::uint32_t offset = static_cast<std::uint32_t>(out.size());
      // local file header
      put32(out, 0x04034b50u);
      header(out, m);
      out += m.name;
      out += m.packed;
      // central directory entry
      put32(central, 0x02014b50u);
      put16(central, 20); // version made by
      header(central, m);
      put16(central, 0); // comment length
      put16(central, 0); // disk number
      put16(central, 0); // internal attributes
      put32(central, 0); // external attributes
      put32(central, offset);
      central += m.name;
    }
    const std::uint32_t cd_offset = static_cast<std::uint32_t>(out.size());
    out += central;
    put32(out, 0x06054b50u); // end of central directory
    put16(out, 0);
    put16(out, 0);
    put16(out, static_cast<std::uint16_t>(m_members.size()));
    put16(out, static_cast<std::uint16_t>(m_members.size()));
    put32(out, static_cast<std::uint32_t>(central.size()));
    put32(out, cd_offset);
    put16(out, 0);

    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (f == nullptr) {
      throw std::runtime_error("NpzWriter: cannot open " + path + " for writing");
    }
    const bool ok = std::fwrite(out.data(), 1, out.size(), f) == out.size();
    if (std::fclose(f) != 0 || !ok) {
      throw std::runtime_error("NpzWriter: error writing " + path);
    }
  }

private:
  struct Member {
    std::string   name;   // with .npy
    std::string   packed; // member bytes as stored in the archive
    std::uint32_t crc    = 0;
    std::uint32_t size   = 0; // uncompressed
    std::uint16_t method = 0; // 0 stored, 8 deflated
  };
  std::vector<Member> m_members;

  template <typename T> static std::string descr()
  {
    const char kind = std::is_floating_point_v<T> ? 'f' : (std::is_signed_v<T> ? 'i' : 'u');
    return std::string(sizeof(T) == 1 ? "|" : "<") + kind + std::to_string(sizeof(T));
  }

  static void put16(std::string& s, std::uint16_t v)
  {
    s.push_back(static_cast<char>(v & 0xff));
    s.push_back(static_cast<char>(v >> 8));
  }
  static void put32(std::string& s, std::uint32_t v)
  {
    for (int i = 0; i < 4; ++i) {
      s.push_back(static_cast<char>((v >> (8 * i)) & 0xff));
    }
  }

  // fields shared by the local header and the central directory entry
  static void header(std::string& s, const Member& m)
  {
    put16(s, 20);           // version needed
    put16(s, 0x0800);       // flags: names are UTF-8
    put16(s, m.method);     // compression
    put16(s, 0);            // time
    put16(s, 0x0021);       // date: 1980-01-01 (the archive is reproducible)
    put32(s, m.crc);
    put32(s, static_cast<std::uint32_t>(m.packed.size()));
    put32(s, m.size);
    put16(s, static_cast<std::uint16_t>(m.name.size()));
    put16(s, 0);            // extra field length
  }

  static std::uint32_t crc32_of(const std::string& data)
  {
#ifdef SOLID_HAVE_ZLIB
    return static_cast<std::uint32_t>(::crc32(0L, reinterpret_cast<const Bytef*>(data.data()), static_cast<uInt>(data.size())));
#else
    static const std::array<std::uint32_t, 256> table = [] {
      std::array<std::uint32_t, 256> t{};
      for (std::uint32_t i = 0; i < 256; ++i) {
        std::uint32_t c = i;
        for (int k = 0; k < 8; ++k) {
          c = (c & 1u) ? 0xedb88320u ^ (c >> 1) : c >> 1;
        }
        t[i] = c;
      }
      return t;
    }();
    std::uint32_t c = 0xffffffffu;
    for (unsigned char b : data) {
      c = table[(c ^ b) & 0xffu] ^ (c >> 8);
    }
    return c ^ 0xffffffffu;
#endif
  }

  static std::u32string decodeUtf8(const std::string& s)
  {
    std::u32string out;
    for (std::size_t i = 0; i < s.size();) {
      const unsigned char c = static_cast<unsigned char>(s[i]);
      int                 n = c < 0x80 ? 0 : (c >> 5) == 0x6 ? 1 : (c >> 4) == 0xe ? 2 : (c >> 3) == 0x1e ? 3 : -1;
      char32_t            u = n == 0 ? c : n == 1 ? (c & 0x1f) : n == 2 ? (c & 0x0f) : (c & 0x07);
      if (n < 0 || i + n >= s.size() + (n == 0 ? 1 : 0)) {
        out.push_back(U'�'); // invalid byte
        ++i;
        continue;
      }
      bool good = true;
      for (int k = 1; k <= n; ++k) {
        const unsigned char d = static_cast<unsigned char>(s[i + k]);
        if ((d >> 6) != 0x2) {
          good = false;
          break;
        }
        u = (u << 6) | (d & 0x3f);
      }
      if (!good) {
        out.push_back(U'�');
        ++i;
        continue;
      }
      out.push_back(u);
      i += n + 1;
    }
    return out;
  }

  void addRaw(const std::string& name, const std::string& descr, const std::vector<std::size_t>& shape, const void* data,
              std::size_t nbytes, bool scalar = false)
  {
    // .npy version 1.0: magic, version, little-endian header length, python dict padded to a multiple of 64 bytes
    std::string dict = "{'descr': '" + descr + "', 'fortran_order': False, 'shape': (";
    if (!scalar) {
      for (std::size_t i = 0; i < shape.size(); ++i) {
        dict += std::to_string(shape[i]) + (shape.size() == 1 ? "," : (i + 1 < shape.size() ? ", " : ""));
      }
    }
    dict += "), }";
    const std::size_t unpadded = 10 + dict.size() + 1;
    dict.append((64 - unpadded % 64) % 64, ' ');
    dict.push_back('\n');
    std::string npy("\x93NUMPY\x01\x00", 8);
    put16(npy, static_cast<std::uint16_t>(dict.size()));
    npy += dict;
    npy.append(static_cast<const char*>(data), nbytes);
    if (npy.size() > 0xfff00000u) {
      throw std::runtime_error("NpzWriter: member " + name + " larger than 4 GB");
    }

    Member m;
    m.name = name + ".npy";
    m.crc  = crc32_of(npy);
    m.size = static_cast<std::uint32_t>(npy.size());
#ifdef SOLID_HAVE_ZLIB
    z_stream zs;
    std::memset(&zs, 0, sizeof(zs));
    if (deflateInit2(&zs, 6, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) == Z_OK) { // raw deflate, as zip wants it
      std::string buf(deflateBound(&zs, static_cast<uLong>(npy.size())), '\0');
      zs.next_in        = reinterpret_cast<Bytef*>(npy.data());
      zs.avail_in       = static_cast<uInt>(npy.size());
      zs.next_out       = reinterpret_cast<Bytef*>(buf.data());
      zs.avail_out      = static_cast<uInt>(buf.size());
      const int status = deflate(&zs, Z_FINISH);
      const auto n     = zs.total_out;
      deflateEnd(&zs);
      if (status == Z_STREAM_END && n < npy.size()) {
        buf.resize(n);
        m.packed = std::move(buf);
        m.method = 8;
      }
    }
#endif
    if (m.method == 0) {
      m.packed = std::move(npy);
    }
    m_members.push_back(std::move(m));
  }
};

} // namespace solid

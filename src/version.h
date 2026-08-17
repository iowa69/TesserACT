// Version and attribution, in one place.
//
// These lived in three: main.cpp for `tesseract-asm --version`, model_main.cpp for
// `tesseract-model --version`, and a literal inside Assembler::run for the version written
// into report.json and report.html. Bumping the first left the other two behind, silently,
// through both 1.2.0 and 1.2.1 -- so a 1.2.1 binary produced reports claiming 1.1.0 and a
// tesseract-model that introduced --marker-density still called itself 1.1.0.
//
// That is worse than cosmetic: the version in a report is provenance. It is the field
// someone reads to find out which code produced a result, and it was wrong by two releases.
#pragma once

namespace ts {

constexpr const char* kVersion = "1.2.4";
constexpr const char* kAuthor = "Giovanni Lorenzin";
constexpr const char* kOrg = "IOWA-BioTech";

}  // namespace ts

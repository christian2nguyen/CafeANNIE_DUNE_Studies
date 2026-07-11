#pragma once

#include "CAFAna/Core/ISyst.h"
#include "CAFAna/Analysis/XSecSystList.h"

#include "StandardRecord/StandardRecord.h"

#include <cassert>
#include <cmath>
#include <iostream>

namespace ana {

enum class XSecDetectorComponent { kBoth, kNear, kFar };

class XSecSyst : public ISyst {
public:
  virtual ~XSecSyst(){};

  void FakeDataDialShift(double sigma, Restorer &restore,
                         caf::StandardRecord *sr, double &weight) const;

  void Shift(double sigma, Restorer &restore, caf::StandardRecord *sr,
             double &weight) const override;

protected:
  XSecSyst(int syst_id, bool applyPenalty = true);
  XSecSyst(int syst_id, const std::string &shortName,
           const std::string &latexName, bool applyPenalty,
           XSecDetectorComponent detComponent);

  friend std::vector<const ISyst *> GetXSecSysts(std::vector<std::string>, bool);
  friend std::vector<const ISyst *>
  GetSplitXSecSysts(std::vector<std::string>, bool);

  int fID;
  XSecDetectorComponent fDetComponent = XSecDetectorComponent::kBoth;
};

std::vector<const ISyst *>
GetXSecSysts(std::vector<std::string> names = {}, bool applyPenalty = true);

std::vector<const ISyst *>
GetSplitXSecSysts(std::vector<std::string> names = {}, bool applyPenalty = true);

} // namespace ana

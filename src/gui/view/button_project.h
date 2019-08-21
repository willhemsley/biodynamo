// Original Author: Bertrand Bellenot   22/08/02
// Modified by: Lukasz Stempniewicz 21/08/19

// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#ifndef GUI_BUTTON_PROJECT_H
#define GUI_BUTTON_PROJECT_H

#include <TGFrame.h>
#include <TGButton.h>

namespace gui {

class ButtonProjectFrame: public TGCompositeFrame {
public:
   // Constructor & destructor
   ButtonProjectFrame(const TGWindow* p, TGWindow* buttonHandler, Int_t newProjectId, Int_t loadProjectId);
   virtual ~ButtonProjectFrame();

   virtual void SetState(int state);

private:
   std::unique_ptr<TGLayoutHints>  fButtonLayout;        // Buttons layout
   std::unique_ptr<TGButton>       fCreateProjectButton; // "Create New Project" button
   std::unique_ptr<TGButton>       fLoadProjectButton;   // "Load Project" button
};

}  // namespace gui

#endif // GUI_BUTTON_PROJECT_H

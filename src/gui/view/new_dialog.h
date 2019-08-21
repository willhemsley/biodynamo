// Author: Lukasz Stempniewicz 21/08/19

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

#ifndef GUI_NEW_DIALOG_H
#define GUI_NEW_DIALOG_H

#include <TColor.h>
#include <TFile.h>
#include <TGFileDialog.h>
#include <TGListBox.h>
#include <TGTab.h>
#include <TGLabel.h>

#include <TGButton.h>
#include <TGFrame.h>
#include <TGTextEntry.h>

#include "core/util/io.h"
#include "core/util/root.h"
#include "gui/view/log.h"
#include "gui/view/model_creator.h"

namespace gui {

extern const char *filetypes[];

class NewProjectDialog : public TGTransientFrame {
 public:
  // Constructor & destructor
  NewProjectDialog(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
                   UInt_t options = kVerticalFrame);
  virtual ~NewProjectDialog();

  // slots
  void              OnCancel();
  void              OnOpen();
  void              CloseWindow();
  Bool_t            OnCreate();
  Bool_t            ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

 private:
  std::unique_ptr<TGCompositeFrame>  fFrame1;
  std::unique_ptr<TGCompositeFrame>  fFrame2;
  std::unique_ptr<TGCompositeFrame>  fFrame3;
  std::unique_ptr<TGVerticalFrame>   fV1;
  std::unique_ptr<TGButton>          fCreateButton;
  std::unique_ptr<TGButton>          fCancelButton;
  std::unique_ptr<TGButton>          fHelpButton;
  std::unique_ptr<TGPictureButton>   fPictButton;
  std::unique_ptr<TGTab>             fTab;
  std::unique_ptr<TGTextEntry>       fTxt1;
  std::unique_ptr<TGTextEntry>       fTxt2;
  std::unique_ptr<TGLabel>           fLerror;
  std::unique_ptr<TGLayoutHints>     fL1;
  std::unique_ptr<TGLayoutHints>     fL2;
  std::unique_ptr<TGLayoutHints>     fL3;
  std::unique_ptr<TGLayoutHints>     fL4;
};

class NewModelDialog : public TGTransientFrame {
 public:
  // Constructor & destructor
  NewModelDialog(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
                 UInt_t options = kVerticalFrame);
  virtual ~NewModelDialog();

  // slots
  void             OnCancel();
  void             CloseWindow();
  Bool_t           OnCreate();
  Bool_t           ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

 private:
  std::unique_ptr<TGCompositeFrame> fFrame1;
  std::unique_ptr<TGCompositeFrame> fFrame2;
  std::unique_ptr<TGVerticalFrame>  fV1;
  std::unique_ptr<TGButton>         fCreateButton;
  std::unique_ptr<TGButton>         fCancelButton;
  std::unique_ptr<TGButton>         fHelpButton;
  std::unique_ptr<TGTab>            fTab;
  std::unique_ptr<TGTextEntry>      fTxt1;
  std::unique_ptr<TGLabel>          fLerror;
  std::unique_ptr<TGLayoutHints>    fL1;
  std::unique_ptr<TGLayoutHints>    fL2;
  std::unique_ptr<TGLayoutHints>    fL3;
  std::unique_ptr<TGLayoutHints>    fL4;
};

}  // namespace gui

#endif // GUI_NEW_DIALOG_H
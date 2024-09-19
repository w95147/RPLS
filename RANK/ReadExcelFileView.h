// ReadExcelFileView.h : interface of the CReadExcelFileView class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_READEXCELFILEVIEW_H__6A650FC5_06A3_486B_967E_4FD34997D7E9__INCLUDED_)
#define AFX_READEXCELFILEVIEW_H__6A650FC5_06A3_486B_967E_4FD34997D7E9__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


class CReadExcelFileView : public CFormView
{
protected: // create from serialization only
	CReadExcelFileView();
	DECLARE_DYNCREATE(CReadExcelFileView)

public:
	//{{AFX_DATA(CReadExcelFileView)
	enum { IDD = IDD_READEXCELFILE_FORM };
	CString	m_Filename;
	CString	m_Rank;
	CString	m_Raw;
	long	m_Colnums;
	//}}AFX_DATA
// Attributes
	_ConnectionPtr pConn;
    _RecordsetPtr pRs;
	CString FileExt;
	BOOL  ifhavefile;
	BOOL  IfpConn;
public:
	CReadExcelFileDoc* GetDocument();

	void ReadFromExcel(CString Biao,long YangBenNum,CString Newbiao);
	void Mydataorder(double* order_data,int* order_order,int ordernum);
	void FenleiNa(CString Newbiao,double* linedata,int Num,CString dataname);
// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CReadExcelFileView)
	public:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	virtual void OnInitialUpdate(); // called first time after construct
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnPrint(CDC* pDC, CPrintInfo* pInfo);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CReadExcelFileView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(CReadExcelFileView)
	afx_msg void OnButton1();
	afx_msg void OnButton2();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in ReadExcelFileView.cpp
inline CReadExcelFileDoc* CReadExcelFileView::GetDocument()
   { return (CReadExcelFileDoc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_READEXCELFILEVIEW_H__6A650FC5_06A3_486B_967E_4FD34997D7E9__INCLUDED_)

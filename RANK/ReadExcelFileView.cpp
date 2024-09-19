// ReadExcelFileView.cpp : implementation of the CReadExcelFileView class
//

#include "stdafx.h"
#include "ReadExcelFile.h"

#include "ReadExcelFileDoc.h"
#include "ReadExcelFileView.h"
#include "odbcinst.h"
#include  "afxdb.h"
#include "comdef.h"
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#include <math.h>
/////////////////////////////////////////////////////////////////////////////
// CReadExcelFileView

IMPLEMENT_DYNCREATE(CReadExcelFileView, CFormView)

BEGIN_MESSAGE_MAP(CReadExcelFileView, CFormView)
	//{{AFX_MSG_MAP(CReadExcelFileView)
	ON_BN_CLICKED(IDC_BUTTON1, OnButton1)
	ON_BN_CLICKED(IDC_BUTTON2, OnButton2)
	//}}AFX_MSG_MAP
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CFormView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, CFormView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CFormView::OnFilePrintPreview)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CReadExcelFileView construction/destruction

CReadExcelFileView::CReadExcelFileView()
	: CFormView(CReadExcelFileView::IDD)
{
	//{{AFX_DATA_INIT(CReadExcelFileView)
	m_Filename = _T("");
	m_Rank = _T("");
	m_Raw = _T("");
	m_Colnums = 0;
	//}}AFX_DATA_INIT
	// TODO: add construction code here
	ifhavefile=FALSE;
	IfpConn=FALSE;
}

CReadExcelFileView::~CReadExcelFileView()
{
}

void CReadExcelFileView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CReadExcelFileView)
	DDX_Text(pDX, IDC_FILENAME, m_Filename);
	DDX_Text(pDX, IDC_RANK, m_Rank);
	DDX_Text(pDX, IDC_RAW, m_Raw);
	DDX_Text(pDX, IDC_COLNUMS, m_Colnums);
	//}}AFX_DATA_MAP
}

BOOL CReadExcelFileView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CFormView::PreCreateWindow(cs);
}

/////////////////////////////////////////////////////////////////////////////
// CReadExcelFileView printing

BOOL CReadExcelFileView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CReadExcelFileView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void CReadExcelFileView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}

void CReadExcelFileView::OnPrint(CDC* pDC, CPrintInfo* /*pInfo*/)
{
	// TODO: add customized printing code here
}

/////////////////////////////////////////////////////////////////////////////
// CReadExcelFileView diagnostics

#ifdef _DEBUG
void CReadExcelFileView::AssertValid() const
{
	CFormView::AssertValid();
}

void CReadExcelFileView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}

CReadExcelFileDoc* CReadExcelFileView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CReadExcelFileDoc)));
	return (CReadExcelFileDoc*)m_pDocument;
}
#endif //_DEBUG

void CReadExcelFileView::OnInitialUpdate()
{
	CFormView::OnInitialUpdate();
	GetParentFrame()->RecalcLayout();
	ResizeParentToFit();
/////////////////////////////////////////////////////////////////////////
	if(FAILED(pConn.CreateInstance(_uuidof(Connection))))
	{
		AfxMessageBox("ADO初始化失败");
		return; 
	} 

}

/////////////////////////////////////////////////////////////////////////////
//ReadFromExcel():读取数据文件
//m_Raw:原数据表名；
//m_Colnums：m_Raw表中的列数 ；
//m_Rank:秩的数据集表名
void CReadExcelFileView::ReadFromExcel(CString Biao,long YangBenNum,CString Newbiao)
{
    variant_t RecordsAffected,vID;
	long totalnum=0;
	CComBSTR  m_bstr;	

	m_bstr.Empty();				
	m_bstr.Append("select count(*) as A0 from [");
	m_bstr.Append(Biao);
	m_bstr.Append("$]");
	try
	 {    
	
		pRs=pConn->Execute((_bstr_t)m_bstr,&RecordsAffected,adCmdText);				   
	}
	catch(_com_error error)
	 {
		CString errorMessage;
		errorMessage.Format("%s",(LPTSTR)error.Description());
		AfxMessageBox(errorMessage);
	 }	
	////////////////////////////////////
	if(!pRs->ADOEOF)
	{
		vID.Clear();
		vID = pRs->GetCollect("A0");
		if(vID.vt!=VT_NULL)	
			totalnum=vID.lVal;		
	}
	pRs->Close();
	if(totalnum==0)
	{
		AfxMessageBox(_T("数据文件中没有有效数据"));
		return;
	}
//////////分配存储数据的空间
	long     datanum=totalnum*YangBenNum;	
	CString*  dataname=new CString[totalnum];
	double*   filedata=new double[datanum];

	for(int k=0;k<datanum;k++)
		filedata[k]=0;
/////////////////////////////////////////
	m_bstr.Empty();				
	m_bstr.Append("select *  from [");
	m_bstr.Append(Biao);
	m_bstr.Append("$]");
	_bstr_t   retu2(m_bstr,FALSE);
	try
	 {    
	
		pRs=pConn->Execute(retu2,&RecordsAffected,adCmdText);				   
	}
	catch(_com_error error)
	 {
		CString errorMessage;
		errorMessage.Format("%s",(LPTSTR)error.Description());
		AfxMessageBox(errorMessage);
	 }	
	////////////////////////////////////
	long  num=0;
	double kr;
	CString tempstr;
	while(!pRs->ADOEOF)
	{
		if(num>=totalnum)
			break;
		vID.Clear();
		vID = pRs->GetCollect(_variant_t(_variant_t((long)0)));//第一列为名字，
		{
			if(vID.vt!=VT_NULL)	
			{
				if(vID.vt==8)
				{
					dataname[num]=vID.bstrVal;
				}
				else
				{
					kr=vID.dblVal;
					tempstr.Format("%.0f",kr);
					dataname[num]=tempstr;
				}		
			}
			else
			{
				pRs->MoveNext();
				continue;
			}
		}
		for(long i=0;i<YangBenNum;i++)
		{
			vID.Clear();
			long  ii=i+1;//.第一列为名字，所以从第2列开始是数据
			vID = pRs->GetCollect(_variant_t(ii));
			if(vID.vt!=VT_NULL)	
				filedata[num*YangBenNum+i]=vID.dblVal;
			else
			{
				filedata[num*YangBenNum+i]=0;
			}
		}
		num++;
		pRs->MoveNext();
	}	
	pRs->Close();
	////////////////////////////////////////////////////////////create new table
	//////////////////////////////////
	m_bstr.Empty();	
	m_bstr.Append("CREATE TABLE [");
	m_bstr.Append(Newbiao);
	m_bstr.Append("] (NAME VARCHAR(50)");

	for(int kl=0;kl<YangBenNum;kl++)
	{
		tempstr.Format(",A%d",kl+1);
		m_bstr.Append(tempstr);
		m_bstr.Append("  int");
	}
	m_bstr.Append(")");
	try
	 {    
	
		pConn->Execute((_bstr_t)m_bstr,&RecordsAffected,adCmdText);				   
	}
	catch(_com_error error)
	 {
		CString errorMessage;
		errorMessage.Format("%s",(LPTSTR)error.Description());
	//	AfxMessageBox(errorMessage);
	 }
	/////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////
	double*  linedata=new double[YangBenNum];
	CString  strW,str1,str2,str3;
	for(int i=0;i<totalnum;i++)
	{
		for(k=0;k<YangBenNum;k++)
		{
			linedata[k]=0;
			linedata[k]=filedata[i*YangBenNum+k];			
		}	
		FenleiNa(Newbiao,linedata,YangBenNum,dataname[i]);
	}
	AfxMessageBox("CONVERSION OK");
}
//CReadExcelFileView：数据排序，并存储文件中
//Newbiao:秩的数据集表的名称
//linedata:行数据
//dataname:行名称
//Num:列总数//
void CReadExcelFileView::FenleiNa(CString Newbiao,double*linedata ,int Num,CString dataname)
{
	variant_t RecordsAffected;	
	CComBSTR  m_bstr;
	long    ALL_num=Num;
	double*  orderdata=new double[ALL_num];
	int*     order=new int [ALL_num];

	for(int k=0;k<ALL_num;k++)
	{
		orderdata[k]=0;
		order[k]=0;
	}
	for(k=0;k<ALL_num;k++)
	{
		orderdata[k]=linedata[k];
	}
///////////数据排序
	Mydataorder(orderdata,order,ALL_num);

	CString  kstr;
	m_bstr.Empty();	
		m_bstr.Append("insert into [");
		m_bstr.Append(Newbiao);
		m_bstr.Append("$] values ('");
	m_bstr.Append(dataname);
	m_bstr.Append("' ");
	for(k=0;k<ALL_num;k++)
	{
		kstr.Format(_T(", %d"),order[k]);
		m_bstr.Append(kstr);
	} 
	m_bstr.Append(")");
	try
	 {    
	
		pConn->Execute((_bstr_t)m_bstr,&RecordsAffected,adCmdText);				   
	}
	catch(_com_error error)
	 {
		CString errorMessage;
		errorMessage.Format("%s",(LPTSTR)error.Description());
		AfxMessageBox(errorMessage);
	 }	
	
}
//数据排序，返回数据的顺序列表
//Mydataorder:需要排序的原数据
//order_order:返回值，行数据对应的秩
//ordernum：数据数目
void CReadExcelFileView::Mydataorder(double* order_data,int* order_order,int ordernum)
{
	double tempshu=0;
	int    tempnum=1;
	int    totalSum=0;
	for(int k=0;k<ordernum;k++)
	{
		tempshu=order_data[k];
		tempnum=1;
		for(int i=0;i<ordernum;i++)
		{
			if(tempshu>order_data[i])
				tempnum++;
		}
		order_order[k]=tempnum;
	}

}
//选择数据文件  
//m_Filename：数据文件完整路径；
//FileExt：文件名扩展。
//本例程只适合处理.xls文件，不能处理.xlsx文件
void CReadExcelFileView::OnButton1() 
{
	CFileDialog dlg(TRUE,_T(".xls"),NULL,OFN_OVERWRITEPROMPT,_T("Excel2003 Files (*.xls)|*.xls|Excel2007 Files (*.xlsx)|*.xlsx||"),NULL);
    if(dlg.DoModal()==IDOK)
    {
        BeginWaitCursor();
        FileExt=dlg.GetFileExt();
        m_Filename=dlg.GetPathName();		        	
	}
	ifhavefile=TRUE;
	UpdateData(FALSE);
}
//CONVERSION
void CReadExcelFileView::OnButton2() 
{
	// TODO: Add your control notification handler code here
	UpdateData(TRUE);
	if(IfpConn)
		pConn->Close();
	if(!ifhavefile)
	{
		AfxMessageBox("Pleas select converted file!");
		return;
	}
	CString   adoinfo; 
    if(FileExt=="xlsx")
        adoinfo.Format(_T("Provider=Microsoft.ACE.OLEDB.12.0;Data Source=%s;Extended Properties=Excel 12.0Xml"),m_Filename);
    else
        adoinfo.Format(_T("Provider=Microsoft.Jet.OLEDB.4.0;Data Source=%s;Extended Properties=Excel 8.0"),m_Filename);
	if(FAILED(pConn->Open(_bstr_t(adoinfo),"","",0)))
	{
		AfxMessageBox("ADO初始化失败");
		return;
	}
	else
		IfpConn=TRUE;
	if(m_Raw.GetLength()==0)
	{
		AfxMessageBox("[raw data sheet]  NO Sheet name!");
		return;
	}
/*	if(m_Rank.GetLength()==0)
	{
		AfxMessageBox("[Rank sheet]  NO Sheet name!");
		return;
	}*/
	if(m_Colnums==0)
	{
		AfxMessageBox("[Sheet  columns]  NO Error!");
		return;
	}
	m_Rank=m_Raw+"Rank";
	ReadFromExcel(m_Raw,m_Colnums-1,m_Rank);
}

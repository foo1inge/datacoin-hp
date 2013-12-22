#include "miningdialog.h"
#include "ui_miningdialog.h"

#include "util.h"
#include "walletmodel.h"

MiningDialog::MiningDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MiningDialog),
    model(NULL)
{
    ui->setupUi(this);
}

MiningDialog::~MiningDialog()
{
    delete ui;
}

void addRowToTable(QTableWidget &table, int n, const char* pstr) {
    QTableWidgetItem *itemRowHeader = table.verticalHeaderItem(n);
    if (NULL == itemRowHeader) itemRowHeader = new QTableWidgetItem;
    itemRowHeader->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
    itemRowHeader->setText(pstr);
    table.setVerticalHeaderItem(n, itemRowHeader);

    table.verticalHeader()->setResizeMode(QHeaderView::Fixed);
}

void MiningDialog::setModel(WalletModel *model)
{
    this->model = model;

    startTimer(1000);

    if (model->getMining())
        ui->startButton->setText("Stop mining");
    else
        ui->startButton->setText("Start mining");

    ui->paramsTable->setRowCount(6);
    ui->paramsTable->setColumnCount(1);

    addRowToTable(*ui->paramsTable, 0, "Mining threads");
    addRowToTable(*ui->paramsTable, 1, "Block count");
    addRowToTable(*ui->paramsTable, 2, "Difficulty");
    addRowToTable(*ui->paramsTable, 3, "Chains per day");
    addRowToTable(*ui->paramsTable, 4, "Chains per minute");
    addRowToTable(*ui->paramsTable, 5, "Primes per second");

    ui->paramsTable->resizeRowsToContents();

    QTableWidgetItem *itemColHeader = ui->paramsTable->horizontalHeaderItem(0);
    if (NULL == itemColHeader) itemColHeader = new QTableWidgetItem;
    itemColHeader->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    itemColHeader->setText("Value");
    ui->paramsTable->setHorizontalHeaderItem(0, itemColHeader);
    ui->paramsTable->horizontalHeader()->setStretchLastSection(true);

    ui->threadsSlider->setMaximum(boost::thread::hardware_concurrency());
}

void MiningDialog::on_startButton_clicked()
{
    if (!model) return;

    if (model->getMining()) {
        ui->threadsSlider->setEnabled(true);
        model->setMining(0);
        ui->startButton->setText("Start mining");
    } else {
        ui->threadsSlider->setEnabled(false);
        model->setMining(ui->threadsSlider->value());
        ui->startButton->setText("Stop mining");
    }
}

extern double dPrimesPerSec;
extern double dChainsPerMinute;
extern double dChainsPerDay;

void setValueToCell(QTableWidget &table, int row, int col, std::string value) {
    QTableWidgetItem *item = table.item(row, col);
    if (NULL == item) item = new QTableWidgetItem;
    item->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    item->setText(value.c_str());
    item->setFlags(item->flags() & ~Qt::ItemIsEditable & ~Qt::ItemIsSelectable);
    table.setItem(row, col, item);
}

void setValueToCell(QTableWidget &table, int row, int col, double value) {
    setValueToCell(table, row, col, dtostr(value));
}

void setValueToCell(QTableWidget &table, int row, int col, int value) {
    setValueToCell(table, row, col, itostr(value));
}

void MiningDialog::timerEvent(QTimerEvent *event) {
    if (!model) return;

    if (model->getMining())
        ui->startButton->setText("Stop mining");
    else
        ui->startButton->setText("Start mining");

    setValueToCell(*ui->paramsTable, 0, 0, ui->threadsSlider->value());
    setValueToCell(*ui->paramsTable, 1, 0, model->getBlockCount());
    setValueToCell(*ui->paramsTable, 2, 0, model->getDifficulty());
    setValueToCell(*ui->paramsTable, 3, 0, dChainsPerDay);
    setValueToCell(*ui->paramsTable, 4, 0, dChainsPerMinute);
    setValueToCell(*ui->paramsTable, 5, 0, dPrimesPerSec);

    //if (NULL != ui->paramsTable->item(1,1))
    //    ui->paramsTable->item(1,1)->setText(itostr(dChainsPerDay).c_str());
}

#ifndef MININGDIALOG_H
#define MININGDIALOG_H

#include <QDialog>

namespace Ui {
class MiningDialog;
}

class WalletModel;

class MiningDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit MiningDialog(QWidget *parent = 0);
    ~MiningDialog();

    void setModel(WalletModel *model);
    
protected:
    void timerEvent(QTimerEvent *event);

private slots:
    void on_startButton_clicked();

private:
    WalletModel *model;
    Ui::MiningDialog *ui;
};

#endif // MININGDIALOG_H

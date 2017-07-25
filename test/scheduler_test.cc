#include "scheduler.h"
#include "unistd.h"
#include "gtest/gtest.h"
#include "io_util.h"

#define ROOTFILE "bdmFile.root"

namespace bdm {
namespace scheduler_test_internal {

class TestSchedulerRestore : public Scheduler<Cell<Soa>> {
public:
  TestSchedulerRestore(const std::string& restore) : Scheduler<Cell<Soa>>("", restore) {}
  void Execute() override {
    execute_calls++;
  }

  unsigned execute_calls = 0;
};

class TestSchedulerBackup : public Scheduler<Cell<Soa>> {
public:
  TestSchedulerBackup(const std::string& backup) : Scheduler<Cell<Soa>>(backup, "") {}

  void Execute() override {
    // sleep
    usleep(350000);
    // backup should be created every second -> every three iterations
    if(execute_calls_ % 3 != 0 || execute_calls_ == 0) {
      EXPECT_FALSE(FileExists(ROOTFILE));
    } else {
      EXPECT_TRUE(FileExists(ROOTFILE));
      remove(ROOTFILE);
    }
    execute_calls_++;
  }
  unsigned execute_calls_ = 0;
};

TEST(SchedulerTest, NoRestoreFile) {
  remove(ROOTFILE);

  // start restore validation
  TestSchedulerRestore scheduler("");
  scheduler.Simulate(100);
  EXPECT_EQ(100u, scheduler.execute_calls);
  EXPECT_EQ(0u, ResourceManager<Cell<Soa>>::Get()->GetCells()->size());

  scheduler.Simulate(100);
  EXPECT_EQ(200u, scheduler.execute_calls);
  EXPECT_EQ(0u, ResourceManager<Cell<Soa>>::Get()->GetCells()->size());

  scheduler.Simulate(100);
  EXPECT_EQ(300u, scheduler.execute_calls);
  EXPECT_EQ(0u, ResourceManager<Cell<Soa>>::Get()->GetCells()->size());
}

TEST(SchedulerTest, Restore) {
  remove(ROOTFILE);

  // create backup that will be restored later on
  auto cells = Cell<>::NewEmptySoa();
  cells.push_back(Cell<>());
  SimulationBackup backup(ROOTFILE, "");
  backup.Backup(&cells, 149);
  EXPECT_EQ(0u, ResourceManager<Cell<Soa>>::Get()->GetCells()->size());

  // start restore validation
  TestSchedulerRestore scheduler(ROOTFILE);
  // 149 simulation steps have already been calculated. Therefore, this call
  // should be ignored
  scheduler.Simulate(100);
  EXPECT_EQ(0u, scheduler.execute_calls);
  EXPECT_EQ(0u, ResourceManager<Cell<Soa>>::Get()->GetCells()->size());

  // Restore should happen within this call
  scheduler.Simulate(100);
  //   only 51 steps should be simulated
  EXPECT_EQ(51u, scheduler.execute_calls);
  EXPECT_EQ(1u, ResourceManager<Cell<Soa>>::Get()->GetCells()->size());

  // add element to see if if restore happens again
  ResourceManager<Cell<Soa>>::Get()->GetCells()->push_back(Cell<>());

  // normal simulation - no restore
  scheduler.Simulate(100);
  EXPECT_EQ(151u, scheduler.execute_calls);
  EXPECT_EQ(2u, ResourceManager<Cell<Soa>>::Get()->GetCells()->size());

  remove(ROOTFILE);
}

TEST(SchedulerTest, Backup) {
  remove(ROOTFILE);

  TestSchedulerBackup scheduler(ROOTFILE);
  Param::kBackupEveryXSeconds = 1;
  // one simulation step takes 350 ms -> backup should be created every three
  // steps
  scheduler.Simulate(7);
  remove(ROOTFILE);
}

}  // namespace scheduler_test_internal
}  // namespace bdm

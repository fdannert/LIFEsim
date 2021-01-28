from lifesim.core.modules import ObservationModule


class BaselineAdjustModule(ObservationModule):
    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def observe(self,
                nstar: int):
        pass
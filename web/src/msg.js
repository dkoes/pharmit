//object for sending messages to a window, but only after we receive an ack
function Message(data, w, dest) {
	var curWindow = w;
	var curDest = dest;
	var curMsg = data;
	var isAcked = 0;
	
	function receiveMessage(event)
	{
		if(event.data == "ack2")
		{
			isAcked = 1;
		}
	}
	
	function check() {
		if(isAcked) {
			curWindow.postMessage(curMsg,curDest);
			curDest ="";
			curMsg = "";
			curWindow = null;
			isAcked = 0;
		}
		else if(curWindow) {
			curWindow.postMessage("ack", curDest);
			setTimeout(check, 250);
		}
	}
	
	window.addEventListener("message", receiveMessage);
	w.postMessage("ack",dest);		
	setTimeout(check, 250);
}


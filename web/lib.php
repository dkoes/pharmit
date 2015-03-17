<?php

//configuration variables
$db_user = "pharmit";
$db_host = "192.168.6.4";
$db_host = "localhost";
$db_name = "pharmit";
$debug = 1;


//subroutines shared by create.php and index.php

function fail($msg)
{
	error_log($msg);
	exit(1);
}

//aunthenticate user, return error message, or empty string on success
function login($user, $pass)
{
	global $db_user, $db_host, $db_name;
	$db = new mysqli($db_host, $db_user, "", $db_name);
	if (mysqli_connect_errno())
		fail('MySQL connect', mysqli_connect_error());

	($stmt = $db->prepare('SELECT password, maxprivatedbs, maxprivateconfs, maxdbs, maxconfs FROM users WHERE email=?')) ||
	fail('Prepare users', $db->error);
	$stmt->bind_param('s', $user) || fail('Bind user', $db->error);
	$stmt->execute();
	$stmt->store_result();
	if($stmt->num_rows > 0) { //have valid username
		$correctpass = "";
		$maxprivatedbs = 0; $maxprivateconfs = 0; $maxdbs = 0; $maxconfs = 0;
		$stmt->bind_result($correctpass, $maxprivatedbs, $maxprivateconfs, $maxdbs, $maxconfs) || fail('Bind pass', $db->error);
		if(!$stmt->fetch() && $db->errno)
			fail('Fetch pass', $db->error);

		if($correctpass == $pass) {
			session_regenerate_id();
			$_SESSION['userid']  = $user;
			$_SESSION['maxprivatedbs'] = $maxprivatedbs;
			$_SESSION['maxprivateconfs'] = $maxprivateconfs;
			$_SESSION['maxconfs'] = $maxconfs;
			session_write_close();
			return "";
		} else {
			return "Invalid password for $user";
		}

	} else {
		return "Invalid user";
	}
	return "";
}

//return name of user
function getname($user)
{
	global $db_user, $db_host, $db_name;
	$db = new mysqli($db_host, $db_user, "", $db_name);
	if (mysqli_connect_errno())
		fail('MySQL connect', mysqli_connect_error());

	($stmt = $db->prepare('SELECT name FROM users WHERE email=?')) ||
	fail('Prepare users', $db->error);
	$stmt->bind_param('s', $user) || fail('Bind user', $db->error);
	$stmt->execute();
	$stmt->store_result();
	if($stmt->num_rows > 0) { //have valid username
		$stmt->bind_result($name) || fail('Bind name', $db->error);
		if(!$stmt->fetch() && $db->errno)
			fail('Fetch name', $db->error);

		return $name;

	}
	return "$user";
}

//return number of libraries in progress and completed
function getcnts($user)
{
	global $db_user, $db_host, $db_name;
	$db = new mysqli($db_host, $db_user, "", $db_name);
	if (mysqli_connect_errno())
		fail('MySQL connect', mysqli_connect_error());

	//in progress count
	($stmt = $db->prepare("SELECT COUNT(*) FROM `databases` WHERE email=? AND status != 'Error' AND status != 'Completed'")) ||
		fail('Prepare db', $db->error);
	$stmt->bind_param('s', $user) || fail('Bind user', $db->error);
	$stmt->execute();
	$stmt->store_result();
	$inprogress = 0;
	if($stmt->num_rows > 0) { //have valid username
		$stmt->bind_result($inprogress) || fail('Bind cnt', $db->error);
		if(!$stmt->fetch() && $db->errno)
			fail('Fetch cnt', $db->error);

	}
	//now completed
	($stmt = $db->prepare("SELECT COUNT(*) FROM `databases` WHERE email=? AND status = 'Completed'")) ||
	fail('Prepare db', $db->error);
	$stmt->bind_param('s', $user) || fail('Bind user', $db->error);
	$stmt->execute();
	$stmt->store_result();
	$completed = 0;
	if($stmt->num_rows > 0) { //have valid username
		$stmt->bind_result($completed) || fail('Bind cnt', $db->error);
		if(!$stmt->fetch() && $db->errno)
			fail('Fetch cnt', $db->error);

	}

	return ["inprogress" => $inprogress,
			"completed" => $completed
	];
}


?>
